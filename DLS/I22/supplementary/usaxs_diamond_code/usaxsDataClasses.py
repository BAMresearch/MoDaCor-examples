#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Used in conjunction with usaxsIO and usaxsPlotter, complete data processing for
usaxs scans collected using the Bonse-Hart drop-instrument on I22 at DLS.

Processing steps are based on dataMerge, developed by Brian Pauw at BAM:
(https://github.com/BAMresearch/dataMerge/tree/main)
"""

from attrs import fields
import numpy as np
import scipy.constants as sp
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit, leastsq
import matplotlib.pyplot as plt
from typing import Any, NoReturn
from collections.abc import Iterable
from statsmodels.stats.weightstats import DescrStatsW

from usaxsPlotter import (
    scatteringDataPlotter,
    foldedDataPlotter,
    mergedDataPlotter,
    saxsDataPlotter,
)
from usaxsToolbox import gaussFunction, saturationChecker, tetrammShaper


# Mixin class for making a dict-like object out of an attrs class
# from: https://github.com/python-attrs/attrs/issues/879
class gimmeItems:
    """Mixin class to make attrs classes quack like a dictionary (well,
    technically a mutable mapping). ONLY use this with attrs classes.

    Provides keys(), values(), and items() methods in order to be
    dict-like in addition to MutableMapping-like. Also provides pop(),
    but it just raises a TypeError :)
    """

    __slots__ = ()  # May as well save on memory?

    def __iter__(self) -> Iterable:
        for ifield in fields(self.__class__):
            yield ifield.name

    def __len__(self) -> int:
        return len(fields(self.__class__))

    def __getitem__(self, k: str) -> Any:
        """
        Adapted from:
        https://github.com/python-attrs/attrs/issues/487#issuecomment-660727537
        """
        try:
            return self.__getattribute__(k)
        except AttributeError as exc:
            raise KeyError(str(exc)) from None

    def __delitem__(self, v: str) -> NoReturn:
        raise TypeError("Cannot delete fields for attrs classes.")

    def __setitem__(self, k: str, v: Any) -> None:
        self.__setattr__(k, v)

    def pop(self, key, default=None) -> NoReturn:
        raise TypeError("Cannot pop fields from attrs classes.")

    def keys(self) -> Iterable:
        return self.__iter__()

    def values(self) -> Iterable:
        for key in self.__iter__():
            yield self.__getattribute__(key)

    def items(self) -> Iterable:
        for key in self.__iter__():
            yield key, self.__getattribute__(key)


# end copy


class configDataObj(gimmeItems):
    """
    Object that holds all necessary parameters for processing 1D USAXS data.

    Default values are given for some parameters (i.e. input units, energy,
    etc.). Default values assume fly scanning of the downstream axis at 18 keV.
    Data path, scan numbers, background scan numbers, and dark current scan
    numbers are specified to retrieve data.
    Paths are given to necessary datasets, as well as alternate paths (e.g.
    entry/path/to/data and entry1/path/to/data). Each requested data type has
    an attribute (i.e. self.diodeNexusPath) with a list, and each requested
    dataset is a list within this (i.e. if you have two diodes, you will have
    two lists within self.diodeNexusPath). For each of the diodes, you may
    provide as many potential paths to the data as you wish - the input function
    will stop after the first successful data pull - whether it's the correct
    data or not. Data from tetramms has an additional attribute to describe the
    channel the data can be found in (i.e. self.diodeChannels).
    """

    def __init__(self):
        self.dataPath = ""
        self.processedSAXSDataPath = self.dataPath + "/processed"
        self.fileNamePrefix = "i22-"
        self.fileNameSuffix = ".nxs"
        self.scanNos: np.array = []
        self.darkCurrentOffset: np.array = [
            [3.02e-9,
            2e-11,],
            [3.02e-9,
            2e-11,],
        ]  # default offsets if no dark curent can be used
        self.backgroundNos: np.array = []
        self.requestedFrames: float = (2000.0, 20000.0)

        self.titleNexusPath = [
            ["entry1/title", "entry/sample/title", "entry/sample/name"]
        ]
        self.yawNexusPath = [["entry/BSDIODES/dfyaw", "entry/instrument/downstream_yaw/value"]]
        self.diodeNexusPath = [
            ["entry1/user_tetramm/data", "entry/HG_USAXS_DIODE/data"],
            ["entry1/bsdiodes/data", "entry/BSDIODES/data"],
        ]
        self.diodeChannels = [3, 3]
        self.diodeNum = len(self.diodeChannels)
        self.diodeScaleFactors = [[0, 1.0], [0.0, 5.33]]
        self.i0NexusPath = [["entry1/I0/data", "entry/I0/data"]]
        self.i0Channels = [6]
        self.countNexusPath = [
            [
                "entry1/instrument/bsdiodes/count_time",
                "entry/instrument/BSDIODES/count_time",
            ]
        ]

        self.inputUnit = "urad"
        self.inputConverter = 1e-6
        self.outputUnit = "1/Angstroms"
        self.outputConverter = 1e10

        self.energy = 18 / 6.2415064799632e15  # convert keV to Joules
        self.fwhm = 8.99  # FWHM of DS rocking curve in urad
        self.qMinTheo = (
            (4 * np.pi * self.energy)
            / (sp.h * sp.c)
            * np.sin(((self.fwhm) * self.inputConverter) / 2)
            / self.outputConverter
        )
        self.qMin = self.qMinTheo - (self.qMinTheo * 0.5)  # Take a bit off
        self.qMax: float = 8.5e-3
        self.nBins: int = 1024
        self.qArray: np.array = []
        self.binEdges: np.array = []

        self.qCalibrationMultipler: float = 1.001169456

        self.usaxsGainScaling: np.array = [0.0, 1.0]
        self.desmearedScaling: np.array = [0.0, 1.0]
        self.scalarMultipler: float = 1e5

        # Temporary variables, to be removed once desmearing is integrated.
        self.desmearedUSAXSDataPath = self.dataPath + "/processing/desmearedData"

        self.feedback = False

    def setDataPath(self, dataPath):
        self.dataPath = dataPath
        self.setSAXSDataPath(dataPath + "/processed")
        self.setDesmearedDataPath(dataPath + "/processing/desmearedData")
        return

    def setSAXSDataPath(self, dataPath):
        self.processedSAXSDataPath = dataPath
        return

    def setDesmearedDataPath(self, dataPath):
        self.desmearedUSAXSDataPath = dataPath
        return

    def setScanNos(self, scanNos):
        self.scanNos = scanNos
        return

    def setDarkCurrentOffset(self, offsets):
        self.darkCurrentOffset = offsets
        return

    def setBackgroundNos(self, backgroundNos):
        self.backgroundNos = backgroundNos
        return

    def setRequestedFrames(self, requestedFrames):
        self.requestedFrames = requestedFrames
        return

    def setUnits(self, inputUnit, outputUnit):
        """
        Set input and output units for data processing. Automatically
        calculates conversion for future use.

        Fly scans currently give urad, while step scans give mrad.
        """
        self.inputUnit = inputUnit
        if inputUnit == "urad":
            self.inputConverter = 1e-06
        elif inputUnit == "microrad":
            self.inputConverter = 1e-06
        elif inputUnit == "mrad":
            self.inputConverter = 1e-03
        elif inputUnit == "millirad":
            self.inputConverter = 1e-03
        else:
            print(
                "!! Input unit not recognised. Please use one of the following:\nurad\nmicrorad\nmrad\nmillirad"
            )

        self.outputUnit = outputUnit
        if outputUnit == "1/Angstroms":
            self.outputConverter = 1e10
        elif outputUnit == "Angstroms":
            self.outputConverter = 1e10
        elif outputUnit == "1/A":
            self.outputConverter = 1e10
        elif outputUnit == "A":
            self.outputConverter = 1e10
        elif outputUnit == "1/nanometers":
            self.outputConverter = 1e09
        elif outputUnit == "nanometers":
            self.outputConverter = 1e09
        elif outputUnit == "1/nm":
            self.outputConverter = 1e09
        elif outputUnit == "nm":
            self.outputConverter = 1e09
        else:
            print(
                "!! Output unit not recognised. Please use one of the following:\n1/Angstroms\nAngstroms\n1/A\nA\n1/nanometers\nnanometers\n1/nm\nnm"
            )

        return

    def setEnergy(self, energy):
        """
        Set energy (provided in keV), and convert to Joules.
        """
        self.energy = energy / 6.2415064799632e15
        return

    def setFWHM(self, fwhm):
        """
        Set FWHM of the downstream crystal rocking curve (in urad).

        Used to calculate theoretical minimum q.
        """
        self.fwhm = fwhm
        return

    def calcQMin(self):
        """
        Calculate qMin from FWHM of downstream rocking curve (theoretical
        minimum achievable q).
        """
        self.qMinTheo = (
            (4 * np.pi * self.energy)
            / (sp.h * sp.c)
            * np.sin(((self.fwhm) * self.inputConverter) / 2)
            / self.outputConverter
        )
        self.qMin = self.qMinTheo - (self.qMinTheo * 0.5)  # Take a bit off
        return

    def setQMin(self, qMin):
        self.qMin = qMin
        return

    def setQMax(self, qMax):
        self.qMax = qMax
        return

    def setNBins(self, nBins):
        self.nBins = nBins
        return

    def calcQArray(self):
        """
        Calculate q array - linearly spaced between qMin and qMax.
        """
        self.qArray = np.linspace(self.qMin, self.qMax, self.nBins)
        return

    def calcBinEdges(self):
        """
        Calculate bin edges - add one bin either side of qMin & qMax to avoid
        issues with NANs when interpolating.
        """
        be = np.linspace((self.qMin), (self.qMax), (self.nBins + 1))
        self.binEdges = np.concatenate(
            (
                np.array((self.qMin - (be[1] - be[0]),)),
                be,
                np.array((self.qMax + (be[-1] - be[-2]),)),
            )
        )
        return

    def setQCalibrationMultipler(self, mulVal):
        self.qCalibrationMultipler = mulVal
        return

    def setUsaxsGainScaling(self, scaling):
        self.usaxsGainScaling = scaling
        return

    def setDesmearedScaling(self, scaling):
        self.desmearedScaling = scaling
        return

    def setScalarMultipler(self, multiplier):
        self.scalarMultipler = multiplier
        return


class scanDataObj(gimmeItems):
    """
    Object that holds all necessary values for an individual scan.

    Raw data is stored in specific variables, so you may re-process from the
    start should you wish
    """

    def __init__(self):
        # Useful knowledge for processing
        self.scanNo: int = []
        self.requestedFrames: np.array = []
        self.present: bool = False
        self.complete: bool = []
        self.beam: bool = []
        self.saturated: bool = []

        # Data...
        self.countTime: np.array = []
        self.I0Raw: np.array = []  # Reserved for raw data
        self.I0: np.array = []
        self.I0STD: np.array = []
        self.I0SEM: np.array = []  # TODO: implement this...
        self.IRaw: np.array = []  # Reserved for raw data
        self.I: np.array = []
        self.ISTD: np.array = []
        self.ISEM: np.array = []  # TODO: implement this...
        self.IMean: np.array = []
        self.Iwt: np.array = []
        self.yawRaw: np.array = []  # Reserved for raw data
        self.yaw: np.array = []
        self.Q: np.array = []
        self.QSTD: np.array = []
        self.QSEM: np.array = []
        self.Mask: np.array = []
        self.binnedN: np.array = []

    def setScanNo(self, scanNo):
        self.scanNo = scanNo
        return

    def setI0Raw(self, I0):
        self.I0Raw = I0
        return

    def setI0(self, I0):
        self.I0 = I0
        return

    def setIRaw(self, IRaw):
        self.IRaw = IRaw
        return

    def setI(self, I):
        self.I = I
        return

    def setYawRaw(self, yawRaw):
        self.yawRaw = yawRaw
        return

    def setYaw(self, yaw):
        self.yaw = yaw
        return

    def setQ(self, q):
        self.q = q
        return

    def checkProcessing(self):
        """
        Check if self.I, self.yaw and self.I0 are populated, if not populate
        by copying from raw data sets self.IRaw, self.yawRaw and self.I0Raw.
        """
        if len(self.I) == 0:
            self.setI(self.IRaw)
        if len(self.yaw) == 0:
            self.setYaw(self.yawRaw)
        if len(self.I0Raw) >= 1:
            if len(self.I0) == 0:
                self.setI0(self.I0Raw)
        return

    def I0Shaper(self):
        """
        Extract the correct tetramm channel and average points for I0.
        """
        self.checkProcessing()
        A, sA, _ = tetrammShaper(tetrammData=self.I0)

        self.I0 = A
        self.I0STD = sA
        return

    def IShaper(self):
        """
        Extract the correct tetramm channel and average points for I.
        """
        self.checkProcessing()
        A, sA, wA = tetrammShaper(tetrammData=self.I)

        self.I = A
        self.ISTD = sA
        self.Iwt = wA
        return

    def yawShaper(self, treatFirstValue = False):
        """
        Extract the downstream channel for yaw.
        """
        # TODO: add in potential for getting upstream values?
        self.checkProcessing()
        yawDims = np.array(self.yaw).ndim
        if yawDims == 1:
            A = self.yaw
        else:
            A = self.yaw
            print("Unexpected shape of yaw.")
        if treatFirstValue:
            A[0] = np.nan
        self.setYaw(A)
        return

    def shapeCheck(self):
        """
        Check data to be processed is 1D, returns False if shape is not 1D.
        """
        shapeChecker = True
        if self.I[0].ndim != 1:
            shapeChecker = False
        if len(self.yaw) >= 1:  # Not necessary for DC scans
            if self.yaw.ndim != 1:
                shapeChecker = False
        if len(self.I0) >= 1:  # Not necessary for DC scans
            if len(self.I0[0]) >= 1:
                if self.I0[0].ndim != 1:
                    shapeChecker = False
        return shapeChecker

    def setIStats(self):
        """
        Set I statistics - used for determining I mean, etc., for dark current
        scans.
        """
        self.checkProcessing()
        if self.shapeCheck():
            self.IMean = [np.mean(IData) for IData in self.I]
            self.ISTD = [np.std(IData) for IData in self.I]
        else:
            print(
                "!! Dimension of some arrays may be higher than 1. Please reduce data using .IShaper()"
            )
        return

    def lenCheck(self):
        """
        Check I, yaw and I0 (if used) are the same length, returns False if
        lengths differ.

        As issues with triggering I0 are still common, only I0 is compared
        to requestedFrames, or the number of frames which should be in the
        scan. If I0 does not have the coorect number of frames, I0 will not
        be used.
        """
        # TODO: fix I0 and then fix this

        self.checkProcessing()
        if self.shapeCheck():
            if len(self.I0[0]) != self.requestedFrames:
                print(f"!! Emptying I0, as it's length ({len(self.I0[0])} was not equal to the number of requested frames ({self.requestedFrames}))")
                self.I0[0] = []  # Empty I0 to ensure it's not used in the future
                yawLen: int = len(self.yaw)
                ILen: int = len(self.I[0])
                if yawLen == ILen:
                    return True
                else:
                    return False
            elif len(self.I0[0]) == self.requestedFrames:
                yawLen: int = len(self.yaw)
                ILen: int = len(self.I[0])
                I0Len: int = len(self.I0[0])
                if (yawLen == ILen) and (yawLen == I0Len):
                    return True
                else:
                    return False
        else:
            print(
                "!! Dimension of some arrays may be higher than 1. Please reduce data using .IShaper(), .I0Shaper() and .yawShaper()"
            )

    def scanChecker(
        self,
        feedback,
        dcMean,
        dcSTD,
        saturationWidth,
        intensityVariation,
        completenessThresh,
    ):
        """
        Check if each scan can be used.

        Checks for completeness of scan (based on requested frame number),
        presence of beam (based on some simple math), and saturation of
        diode. Returns true if at least one scan can be used for processing
        and false if no scans can be processed.
        """
        self.beam = []
        self.complete = []
        self.saturated = []
        for inc in range(len(self.I)):
            # Assume it won't be usable...
            beam = False
            complete = False
            saturated = False

            if len(self.I[inc]) >= (self.requestedFrames * completenessThresh):
                complete = True
                if feedback:
                    print("Scan holds > 75% of frames! Scan can be used...")
                if max(self.I[inc]) >= (
                    10 * ((3 * dcSTD[inc]) + dcMean[inc])
                ):  # An order of magnitude > 3*stdev of DC + mean of DC
                    beam = True
                    if feedback:
                        print("Beam has been found! Scan can be used...")
                    saturationCheck = saturationChecker(
                        self.I[inc], saturationWidth, intensityVariation
                    )
                    if saturationCheck:
                        if feedback:
                            print(
                                "Diode has been saturated, data can processed but cannot be used for centering."
                            )
                        saturated = True
                    elif not saturationCheck:
                        if feedback:
                            print(
                                "Diode has not been saturated, data can be processed and used for centering."
                            )
                        saturated = False
                elif max(self.I[inc]) < (
                    (10 * np.std(self.I[inc])) + np.mean(self.I[inc])
                ):
                    beam = False
                    if feedback:
                        print("!! No beam was found in scan, cannot use.")

            elif len(self.I[inc]) < (self.requestedFrames * completenessThresh):
                complete = False
                if feedback:
                    print("!! Scan was not complete enough for use.")

            self.beam.append(beam)
            self.complete.append(complete)
            self.saturated.append(saturated)
        return

    def dataTrimmer(self):
        """
        Trim data in I, yaw and I0 (if used) to be the same (smallest) length.

        Only trims I0 if it is the correct length - if tetramm was in free run
        for some reason you don't want to use the data...
        """
        self.checkProcessing()
        if self.shapeCheck():
            print("Shape was okay during trimming")
            if not self.lenCheck():
                yawLen: int = len(self.yaw)
                ILen: int = min([len(iData) for iData in self.I])
                if len(self.I0[0]) != self.requestedFrames:
                    minLen = min((yawLen, ILen))
                elif len(self.I0[0]) == self.requestedFrames:
                    I0Len = len(self.I0[0])
                    minLen = min((yawLen, ILen, I0Len))
                self.yaw = self.yaw[:minLen]
                self.I = [iData[:minLen] for iData in self.I]
                self.ISTD = [iData[:minLen] for iData in self.ISTD]
                self.Iwt = [iData[:minLen] for iData in self.Iwt]
                if len(self.I0[0]) == self.requestedFrames:
                    self.I0[0] = self.I0[0][:minLen]
                    self.I0STD[0] = self.I0STD[0][:minLen]
            else:
                print("!! Length of some arrays was not what was expected.")
            return
        else:
            print(
                "!! Dimension of some arrays may be higher than 1. Please reduce data using .IShaper() or .yawShaper()"
            )

    def timeNorm(self, normI0=True):
        """
        Normalise intensity to count time for each frame.

        Now with the option to ignore I0 during normalisation
        (used for dark current scans).
        """
        for inc in range(len(self.I)):
            if self.countTime.size != len(self.I[inc]):
                countTime = np.full_like(self.I[inc], fill_value=self.countTime)
            else:
                countTime = self.countTime

            self.I[inc] /= countTime
            self.ISTD[inc] /= abs(countTime)

        if normI0:
            if len(self.I0[0]) == self.requestedFrames:
                self.I0[0] /= countTime
                self.I0STD[0] /= abs(countTime)
        return

    def I0Norm(self):
        """
        Normalise intensity to incident intensity for each frame.
        """
        if all(
            (
                len(self.I0[0]) == self.requestedFrames,
                len(self.I[0]) == self.requestedFrames,
            )
        ):
            zeroMask = np.array(self.I0[0] == 0)
            I0 = self.I0[0][~zeroMask]

            for inc in range(len(self.I)):
                IData = self.I[inc][~zeroMask]
                self.I[inc] = IData / I0
                if len(self.ISTD[inc]) >= 1:
                    ISTD = self.ISTD[inc][~zeroMask]
                    self.ISTD[inc] = ISTD / abs(I0)
                self.Iwt[inc] = self.Iwt[inc][~zeroMask]

            self.I0[0] = I0
            self.yaw = self.yaw[~zeroMask]
        else:
            print(
                "!! I0 or I did not record the correct number of frames, no normalisation can be done..."
            )
        return

    def dataSorter(self, side):
        """
        'Fold' and sort data by yaw values.

        To get a sensible return from this function, data should have been
        centred beforehand. Does not sort I0, so if normalisation is needed
        it must be done before sorting.
        """
        self.checkProcessing()
        if self.shapeCheck():
            if self.lenCheck():
                if side == "both":
                    yaw = self.yaw
                elif side == "positive":
                    yaw = self.yaw[self.yaw >= 0]
                elif side == "negative":
                    yaw = self.yaw[self.yaw <= 0]
                else:
                    print(
                        "!! Side of USAXS scan to keep for further processing is not recognised. Please use either:\n'both',\n'positive',\nor 'negative'."
                    )
                inds = abs(yaw.argsort())
                self.yaw = abs(yaw[inds])

                for inc in range(len(self.I)):
                    if side == "both":
                        I = self.I[inc]
                        ISTD = self.ISTD[inc]
                        Iwt = self.Iwt[inc]
                        ISTD = self.ISTD[inc]
                    elif side == "positive":
                        I = self.I[inc][self.yaw > 0]
                        ISTD = self.ISTD[inc][self.yaw > 0]
                        Iwt = self.Iwt[inc][self.yaw > 0]
                    elif side == "negative":
                        I = self.I[inc][self.yaw < 0]
                        ISTD = self.ISTD[inc][self.yaw < 0]
                        Iwt = self.Iwt[inc][self.yaw < 0]

                    self.I[inc] = I[inds]
                    self.ISTD[inc] = ISTD[inds]
                    self.Iwt[inc] = Iwt[inds]
            else:
                print("Data sets are different lengths.")
        else:
            print("!! At least one data set is an unexpected shape.")

    def qConverter(self, configDataObj):
        """
        Convert yaw values to q. Output unit dictated by configDataObj.
        """
        self.Q = (
            configDataObj.qCalibrationMultipler
            * (4 * np.pi * configDataObj.energy)
            / (sp.h * sp.c)
            * np.sin((self.yaw * configDataObj.inputConverter) / 2)
            / configDataObj.outputConverter
        )
        return

class rawUsaxsDataObj(gimmeItems):
    """
    Object that carries necessary data for one run of USAXS data.
    """

    def __init__(self):
        self.runNo: int = []

        self.IT: float = []
        self.peakPos: float = []

        self.loGainDC = scanDataObj()
        self.loGain = scanDataObj()
        self.hiGainDC = scanDataObj()
        self.hiGain = scanDataObj()

    def setRunNo(self, runNo):
        self.runNo = runNo
        return

    def setIT(self, IT):
        self.IT = IT
        return


class processedDataObj(gimmeItems):
    """
    Object that carries necessary data for a processed sample.
    """

    def __init__(self):
        self.dataScanNos: np.array = []
        self.bgScanNos: np.array = []
        self.yaw: np.array = []
        self.yawSTD: np.array = []
        self.yawSEM: np.array = []
        self.I: np.array = []
        self.ISTD: np.array = []
        self.ISEM: np.array = []
        self.IN: np.array = []
        self.Q: np.array = []
        self.QSTD: np.array = []
        self.QSEM: np.array = []
        self.binnedN: np.array = []

        self.desmearedI: np.array = []
        self.desmearedIerr: np.array = []
        self.desmearedQ: np.array = []

        self.saxsScanNo: np.array = []
        self.saxsI: np.array = []
        self.saxsIerr: np.array = []
        self.saxsQ: np.array = []

    def setDataScanNos(self, dataScanNos):
        self.dataScanNos = dataScanNos
        return

    def setBGScanNos(self, bgScanNos):
        self.bgScanNos = bgScanNos
        return

    def setI(self, I, ISTD, ISEM):
        self.I = I
        self.ISTD = ISTD
        self.ISEM = ISEM
        return

    def setYaw(self, yaw, yawSTD, yawSEM):
        self.yaw = yaw
        self.yawSTD = yawSTD
        self.yawSEM = yawSEM
        return

    def setQ(self, Q, QSTD, QSEM):
        self.Q = Q
        self.QSTD = QSTD
        self.QSEM = QSEM
        return


class scatteringDataObj(gimmeItems):
    """
    Object that carries necessary data and processing steps for a sample.

    Includes sample name, low gain, high gain, and merged data. Processing
    steps include shaping of diode data sets (I0, bs3diode, panda), setting
    intensity statistics (for dark current scans), trimming datasets (useful
    for prematurely stops of the mapping queue), dark current subtraction,
    data masking (based on maximum intensity of high gain scan), data centring,
    data sorting (including 'folding' and IT calculation), q conversion,
    binning, scaling, merging, and background subtraction.
    """

    def __init__(self):
        self.sampleName: str = ""
        self.scanNos: list = ()
        self.usaxsRuns: dict = {}
        self.mData = processedDataObj()
        self.processingStep: int = (
            0  # Purely for naming figures when you're plotting...
        )
        self.processingSteps: list = []  # For information purposes, will be inserted into metadata of processed and merged data
        self.swaxsPresent: bool = False
        self.usaxsPresent: bool = False
        self.sharedProcessingDone: bool = False

    def setSampleName(self, sampleName):
        self.sampleName = sampleName
        return

    def setScanNos(self):
        self.scanNos = list(self.scanNos)  # TODO: init this properly in __init__
        for runNo in list(self.usaxsRuns.keys()):
            if self.usaxsRuns[runNo].loGain.present:
                self.scanNos.append(str(self.usaxsRuns[runNo].loGain.scanNo))
            if self.usaxsRuns[runNo].hiGain.present:
                self.scanNos.append(str(self.usaxsRuns[runNo].hiGain.scanNo))
        return

    def setRequestedFrameNos(self, configDataObj):
        for runNo in list(self.usaxsRuns.keys()):
            self.usaxsRuns[
                runNo
            ].loGainDC.requestedFrames = configDataObj.requestedFrames[0]
            self.usaxsRuns[
                runNo
            ].hiGainDC.requestedFrames = configDataObj.requestedFrames[0]
            self.usaxsRuns[
                runNo
            ].loGain.requestedFrames = configDataObj.requestedFrames[1]
            self.usaxsRuns[
                runNo
            ].hiGain.requestedFrames = configDataObj.requestedFrames[1]
        return

    def setScanDetails(self, configDataObj):
        self.setScanNos()
        self.setRequestedFrameNos(configDataObj)
        return

    def setIStats(self):
        """
        Calculate I statistical values for high and low gain scans. Used for dark
        current scans.

        Also determines how many diodes you're using for data collection.
        """
        for runNo in list(self.usaxsRuns.keys()):
            if self.usaxsRuns[runNo].loGainDC.present:
                self.usaxsRuns[runNo].loGainDC.setIStats()
            if self.usaxsRuns[runNo].loGain.present:
                loIN = len(self.usaxsRuns[runNo].loGain.I)
            else:
                loIN = 0
            if self.usaxsRuns[runNo].hiGainDC.present:
                self.usaxsRuns[runNo].hiGainDC.setIStats()
            if self.usaxsRuns[runNo].hiGain.present:
                hiIN = len(self.usaxsRuns[runNo].hiGain.I)
            else:
                hiIN = 0
        self.IN = max(loIN, hiIN)
        return

    def scanChecker(self, configDataObj, **kwargs):
        """
        Checks for validity of scans within the scatteringDataObj.

        Checks high and low gain scans (if present) for completeness (scans are
        considered complete if > 75% of requested frames are present) and
        presence of a peak (rocking curve). If scans are determined to be
        saturated or unavailable, processing will fail.
        """
        saturationWidth = kwargs.get("saturationWidth", 15)
        intensityVariation = kwargs.get("intensityVariation", 0.01)
        completenessThresh = kwargs.get("completenessThresh", 0.75)

        feedback = configDataObj.feedback
        usableScans = 0

        for runNo in list(self.usaxsRuns.keys()):
            if self.usaxsRuns[runNo].loGain.present:
                if self.usaxsRuns[runNo].loGainDC.present:
                    if feedback:
                        print("Checking loGain scan for beam using dark current data")
                    self.usaxsRuns[runNo].loGain.scanChecker(
                        feedback=feedback,
                        dcMean=self.usaxsRuns[runNo].loGainDC.IMean,
                        dcSTD=self.usaxsRuns[runNo].loGainDC.ISTD,
                        saturationWidth=saturationWidth,
                        intensityVariation=intensityVariation,
                        completenessThresh=completenessThresh,
                    )
                else:
                    # TODO: check the use of dark current default values sifts out scans where there is no beam
                    if feedback:
                        print(
                            "Checking loGain scan for beam using default dark current values"
                        )
                    self.usaxsRuns[runNo].loGain.scanChecker(
                        feedback=feedback,
                        dcMean=configDataObj.darkCurrentOffset[0],
                        dcSTD=configDataObj.darkCurrentOffset[0],
                        saturationWidth=saturationWidth,
                        intensityVariation=intensityVariation,
                        completenessThresh=completenessThresh,
                    )
            else:
                self.usaxsRuns[runNo].loGain.beam = []
                self.usaxsRuns[runNo].loGain.complete = []
                self.usaxsRuns[runNo].loGain.saturated = []
                for inc in range(self.IN):
                    self.usaxsRuns[runNo].loGain.beam.append(False)
                    self.usaxsRuns[runNo].loGain.complete.append(False)

            if self.usaxsRuns[runNo].hiGain.present:
                if self.usaxsRuns[runNo].hiGainDC.present:
                    if feedback:
                        print("Checking hiGain scan for beam using dark current data")
                    self.usaxsRuns[runNo].hiGain.scanChecker(
                        feedback,
                        self.usaxsRuns[runNo].hiGainDC.IMean,
                        self.usaxsRuns[runNo].hiGainDC.ISTD,
                        saturationWidth,
                        intensityVariation,
                        completenessThresh,
                    )
                else:
                    # TODO: check the use of dark current default values sifts out scans where there is no beam
                    if feedback:
                        print(
                            "Checking hiGain scan for beam using default dark current values"
                        )
                    self.usaxsRuns[runNo].hiGain.scanChecker(
                        feedback,
                        configDataObj.darkCurrentOffset[1],
                        configDataObj.darkCurrentOffset[1],
                        saturationWidth,
                        intensityVariation,
                        completenessThresh,
                    )
            else:
                self.usaxsRuns[runNo].hiGain.beam = []
                self.usaxsRuns[runNo].hiGain.complete = []
                self.usaxsRuns[runNo].hiGain.saturated = []
                for inc in range(self.IN):
                    self.usaxsRuns[runNo].hiGain.beam.append(False)
                    self.usaxsRuns[runNo].hiGain.complete.append(False)

            if all((self.usaxsRuns[runNo].loGain.beam)):
                usableScans += 1
            if all((self.usaxsRuns[runNo].hiGain.beam)):
                usableScans += 1

        if usableScans >= 1:
            return True
        elif usableScans == 0:
            return False

    def I0Shaper(self, feedback):
        """
        Extract the correct tetramm channel and average points for I0 for low
        gain and high gain scans.

        If I0 has not collected the requested number of frames, data will not
        be processed and no further normalisation to I0 can be done.
        """
        for runNo in list(self.usaxsRuns.keys()):
            if self.usaxsRuns[runNo].loGainDC.present:
                self.usaxsRuns[runNo].loGainDC.I0Shaper()
            else:
                if feedback:
                    print("!! loGainDC I0 was skipped")
            if self.usaxsRuns[runNo].loGain.present:
                self.usaxsRuns[runNo].loGain.I0Shaper()
                if feedback:
                    print(f"loGain I0 was shaped, it was {len(self.usaxsRuns[runNo].loGain.I0[0])} long")
            else:
                if feedback:
                    print("!! loGain I0 was skipped")
            if self.usaxsRuns[runNo].hiGainDC.present:
                self.usaxsRuns[runNo].hiGainDC.I0Shaper()
            else:
                if feedback:
                    print("!! hiGainDC I0 was skipped")
            if self.usaxsRuns[runNo].hiGain.present:
                self.usaxsRuns[runNo].hiGain.I0Shaper()
                if feedback:
                    print(f"hiGain I0 was shaped, it was {len(self.usaxsRuns[runNo].hiGain.I0[0])} long")
            else:
                if feedback:
                    print("!! hiGain I0 was skipped")
        return

    def IShaper(self, feedback):
        """
        Extract the correct tetramm channel and average points for I for low
        and high gain scans.
        """
        for runNo in list(self.usaxsRuns.keys()):
            if self.usaxsRuns[runNo].loGainDC.present:
                self.usaxsRuns[runNo].loGainDC.IShaper()
                if feedback:
                    print(f"loGainDC I was shaped, it was {len(self.usaxsRuns[runNo].loGainDC.I)} long")
            if self.usaxsRuns[runNo].loGain.present:
                self.usaxsRuns[runNo].loGain.IShaper()
                if feedback:
                    print(f"loGain I was shaped, it was {len(self.usaxsRuns[runNo].loGain.I)} long")
            if self.usaxsRuns[runNo].hiGainDC.present:
                self.usaxsRuns[runNo].hiGainDC.IShaper()
                if feedback:
                    print(f"hiGainDC I was shaped, it was {len(self.usaxsRuns[runNo].hiGainDC.I)} long")
            if self.usaxsRuns[runNo].hiGain.present:
                self.usaxsRuns[runNo].hiGain.IShaper()
                if feedback:
                    print(f"hiGain I was shaped, it was {len(self.usaxsRuns[runNo].hiGain.I)} long")
        return

    def yawShaper(self, treatFirstValue, feedback):
        """
        Extract the downstream channel for yaw for low and high gain scans.
        """
        for runNo in list(self.usaxsRuns.keys()):
            if self.usaxsRuns[runNo].loGain.present:
                self.usaxsRuns[runNo].loGain.yawShaper(treatFirstValue=treatFirstValue)
                if feedback:
                    print("loGain yaw was shaped")
            else:
                if feedback:
                    print("!! loGain yaw was skipped")
            if self.usaxsRuns[runNo].hiGain.present:
                self.usaxsRuns[runNo].hiGain.yawShaper(treatFirstValue=treatFirstValue)
                if feedback:
                    print("hiGain yaw was shaped")
            else:
                if feedback:
                    print("!! hiGain yaw was skipped")
        return

    def dataTrimmer(self):
        """
        Trim data in I, yaw and I0 (if used) to be the same (smallest) length,
        for low and high gain scans.

        Data for each scan is trimmed spearately, i.e. if the low gain scan is
        shorter than the high gain scan, the scans will be different lengths.
        Dark current scans are not currently trimmed.
        """
        for runNo in list(self.usaxsRuns.keys()):
            if self.usaxsRuns[runNo].loGain.present:
                self.usaxsRuns[runNo].loGain.dataTrimmer()
            if self.usaxsRuns[runNo].hiGain.present:
                self.usaxsRuns[runNo].hiGain.dataTrimmer()
        return

    def dataShaper(self, **kwargs):
        """
        Fully shapes all aspects of data (i.e. I, I0 and yaw). Plots data if
        requested.

        Shapes I, yaw, I0, and subsequently trims data to ensure all arrays for
        the same scan are 1d and the same length. Channel can be specified for
        I, I0, and yaw. The default for all channels is I = 2 (front diode),
        I0 = 6 (QBPM2_total), and yaw = 2 (downstream). This channel selection
        will only work if the full dataset (i.e. full tetramm or both encoder
        channels on the panda) has been acquired.
        If configDataObj is passed in and plotMe is "True", data is plotted
        and plot is saved to the "/processing/Figures/" subdirectory in the
        specified dataPath.
        """
        treatFirstValue = kwargs.get("treatFirstValue", False)
        configDataObj = kwargs.get("configDataObj", None)
        plotMe = kwargs.get("plotMe", False)
        errBars = kwargs.get("errBars", False)

        if configDataObj is not None:
            feedback = configDataObj.feedback
        else:
            feedback = False

        print("Shaping all data...")

        if len(self.scanNos) >= 1:
            if feedback:
                print("Shaping diode data...")
            self.IShaper(feedback)
            if feedback:
                print("Shaping panda data...")
            self.yawShaper(treatFirstValue=treatFirstValue,feedback=feedback)
            if feedback:
                print("Shaping I0 data...")
            self.I0Shaper(feedback)
            if feedback:
                print("Trimming data...")
            self.dataTrimmer()

            self.setIStats()
            self.processingSteps.append("Shape data (I, I0, yaw)")
        elif len(self.scanNos) == 0:
            print("!! No USAXS scans are present for this sample, perhaps check if you can lazy load this data?")

        # TODO: Fix the naming convention for these plots...
        if plotMe:
            if configDataObj is not None:
                scanNosString = "i22"
                for scanNo in self.scanNos:
                    scanNosString += "-"
                    scanNosString += str(scanNo)
                fileName = (
                    configDataObj.dataPath
                    + "/processing/Figures/"
                    + scanNosString
                    + "_step"
                    + str(self.processingStep)
                    + "_Shaping.png"
                )
                self.processingStep += 1
                scatteringDataPlotter(
                    scatteringDataObject=self,
                    configDataObj=configDataObj,
                    xUnit="yaw",
                    errBars=errBars,
                    mkrSize=3,
                    fileName=fileName,
                    yScale="log",
                    plotTitle=(
                        "Shaped Data: " + scanNosString
                    ),
                )
            else:
                print(
                    "!! Please provide configDataObj if you would like to view the plot."
                )
        return

    def timeNorm(self, **kwargs):
        """
        Normalise data to count time.

        Normalises intensity (I and I0) data to count time found in nexus file.
        If count time was not found in original nexus file, it is set to a
        default of 1.
        If configDataObj is passed in and plotMe is "True", data is plotted
        and plot is saved to the "/processing/Figures/" subdirectory in the
        specified dataPath.
        """
        configDataObj = kwargs.get("configDataObj", None)
        plotMe = kwargs.get("plotMe", False)
        errBars = kwargs.get("errBars", False)

        if configDataObj is not None:
            feedback = configDataObj.feedback
        else:
            feedback = False

        print("Normalising to time...")
        for runNo in list(self.usaxsRuns.keys()):
            if all(
                (
                    self.usaxsRuns[runNo].loGain.complete
                    + self.usaxsRuns[runNo].loGain.beam
                )
            ):
                if self.usaxsRuns[runNo].loGain.present:
                    if feedback:
                        print("Normalising low gain usaxs scan to count time")
                    self.usaxsRuns[runNo].loGain.timeNorm(normI0=True)
                if self.usaxsRuns[runNo].loGainDC.present:
                    if feedback:
                        print("Normalising low gain dark current scan to count time")
                    self.usaxsRuns[runNo].loGainDC.timeNorm(normI0=False)
            else:
                if feedback:
                    print(
                        "!! Low gain scan either incomplete or did not have beam, cannot be normalised to time..."
                    )

            if all(
                (
                    self.usaxsRuns[runNo].hiGain.complete
                    + self.usaxsRuns[runNo].hiGain.beam
                )
            ):
                if self.usaxsRuns[runNo].hiGain.present:
                    if feedback:
                        print("Normalising high gain usaxs scan to count time")
                    self.usaxsRuns[runNo].hiGain.timeNorm(normI0=True)
                if self.usaxsRuns[runNo].hiGainDC.present:
                    if feedback:
                        print("Normalising high gain dark current scan to count time")
                    self.usaxsRuns[runNo].hiGainDC.timeNorm(normI0=False)
            else:
                if feedback:
                    print(
                        "!! High gain scan either incomplete or did not have beam, cannot be normalised to time..."
                    )

        self.setIStats()
        self.processingSteps.append("Normalise data to time")

        if plotMe:
            if configDataObj is not None:
                scanNosString = "i22"
                for scanNo in self.scanNos:
                    scanNosString += "-"
                    scanNosString += str(scanNo)
                fileName = (
                    configDataObj.dataPath
                    + "/processing/Figures/"
                    + scanNosString
                    + "_step"
                    + str(self.processingStep)
                    + "_TimeNorm.png"
                )
                self.processingStep += 1
                scatteringDataPlotter(
                    scatteringDataObject=self,
                    configDataObj=configDataObj,
                    xUnit="yaw",
                    errBars=errBars,
                    mkrSize=3,
                    fileName=fileName,
                    yScale="log",
                    plotTitle=(
                        "Time Normalised Data: " + scanNosString
                    ),
                )
            else:
                print(
                    "!! Please provide configDataObj if you would like to view the plot."
                )
        return

    def darkCurrentSubtraction(self, **kwargs):
        """
        Subtract specified dark current (mean) from high and low gain scans.

        Dark current can be subtracted by one of two methods:
        1) basic - which subtracts a fixed value saved in the configDataObj
        2) interleaved - which subtracts the mean of the interleaved dark
           current scan (collected directly prior to the usaxs scan)

        If configDataObj is passed in and plotMe is "True", data is plotted
        and plot is saved to the "/processing/Figures/" subdirectory in the
        specified dataPath.
        """
        configDataObj = kwargs.get("configDataObj", None)
        plotMe = kwargs.get("plotMe", False)
        errBars = kwargs.get("errBars", False)
        subtractionMethod = kwargs.get("subtractionMethod", "interleaved")

        if configDataObj is not None:
            feedback = configDataObj.feedback
        else:
            feedback = False

        if subtractionMethod == "basic":
            if configDataObj is not None:
                print("Doing a basic dark current subtraction...")
                for runNo in list(self.usaxsRuns.keys()):
                    if all(
                        (
                            self.usaxsRuns[runNo].loGain.complete,
                            +self.usaxsRuns[runNo].loGain.beam,
                        )
                    ):
                        for inc in range(self.IN):
                            A = self.usaxsRuns[runNo].loGain.I[inc]
                            B = configDataObj.darkCurrentOffset[inc][0]
                            self.usaxsRuns[runNo].loGain.I[inc] = A - B
                    else:
                        if feedback:
                            print(
                                "!! Low gain scan either incomplete or did not have beam, dark current subtraction can not be performed..."
                            )

                    if all(
                        (
                            self.usaxsRuns[runNo].hiGain.complete
                            + self.usaxsRuns[runNo].hiGain.beam
                        )
                    ):
                        for inc in range(self.IN):
                            A = self.usaxsRuns[runNo].hiGain.I[inc]
                            B = configDataObj.darkCurrentOffset[inc][1]
                            self.usaxsRuns[runNo].hiGain.I[inc] = A - B
                    else:
                        if feedback:
                            print(
                                "!! High gain scan either incomplete or did not have beam, dark current subtraction can not be performed..."
                            )
            else:
                print(
                    "!! Please provide a configDataObj to access basic dark current subtraction..."
                )
            self.processingSteps.append("Dark current subtraction - basic")

        elif subtractionMethod == "interleaved":
            try:
                print("Doing an interleaved dark current subtraction...")
                for runNo in list(self.usaxsRuns.keys()):
                    if all(
                        (
                            self.usaxsRuns[runNo].loGain.complete
                            + self.usaxsRuns[runNo].loGain.beam
                        )
                    ):
                        for inc in range(self.IN):
                            A = self.usaxsRuns[runNo].loGain.I[inc]
                            B = self.usaxsRuns[runNo].loGainDC.IMean[inc]
                            self.usaxsRuns[runNo].loGain.I[inc] = A - B
                    else:
                        if feedback:
                            print(
                                "!! Low gain scan either incomplete or did not have beam, dark current subtraction can not be performed..."
                            )

                    if all(
                        (
                            self.usaxsRuns[runNo].hiGain.complete
                            + self.usaxsRuns[runNo].hiGain.beam
                        )
                    ):
                        for inc in range(self.IN):
                            A = self.usaxsRuns[runNo].hiGain.I[inc]
                            B = self.usaxsRuns[runNo].hiGainDC.IMean[inc]
                            self.usaxsRuns[runNo].hiGain.I[inc] = A - B
                    else:
                        if feedback:
                            print(
                                "!! High gain scan either incomplete or did not have beam, dark current subtraction can not be performed..."
                            )
                self.processingSteps.append("Dark current subtraction - interleaved")
            except (AttributeError, IndexError):
                if configDataObj is not None:
                    print(
                        "!! An error was encountered doing an interleaved dark current subtraction, now doing a basic dark current subtraction..."
                    )
                    for runNo in list(self.usaxsRuns.keys()):
                        if all(
                            (
                                self.usaxsRuns[runNo].loGain.complete
                                + self.usaxsRuns[runNo].loGain.beam
                            )
                        ):
                            for inc in range(self.IN):
                                A = self.usaxsRuns[runNo].loGain.I[inc]
                                B = configDataObj.darkCurrentOffset[inc][0]
                                self.usaxsRuns[runNo].loGain.I[inc] = A - B
                        else:
                            if feedback:
                                print(
                                    "!! Low gain scan either incomplete or did not have beam, dark current subtraction can not be performed..."
                                )

                        if all(
                            (
                                self.usaxsRuns[runNo].hiGain.complete
                                + self.usaxsRuns[runNo].hiGain.beam
                            )
                        ):
                            for inc in range(self.IN):
                                A = self.usaxsRuns[runNo].hiGain.I[inc]
                                B = configDataObj.darkCurrentOffset[inc][1]
                                self.usaxsRuns[runNo].hiGain.I[inc] = A - B
                        else:
                            if feedback:
                                print(
                                    "!! High gain scan either incomplete or did not have beam, dark current subtraction can not be performed..."
                                )
                else:
                    print(
                        "!! Please provide a configDataObj to access basic dark current subtraction..."
                    )
                self.processingSteps.append("Dark current subtraction - basic")

        else:
            print(
                "!! Dark current subtraction method was not recognised, please specify either basic or interleaved."
            )

        if plotMe:
            if configDataObj is not None:
                scanNosString = "i22"
                for scanNo in self.scanNos:
                    scanNosString += "-"
                    scanNosString += str(scanNo)
                fileName = (
                    configDataObj.dataPath
                    + "/processing/Figures/"
                    + scanNosString
                    + "_step"
                    + str(self.processingStep)
                    + "_DCSub.png"
                )
                self.processingStep += 1
                scatteringDataPlotter(
                    scatteringDataObject=self,
                    configDataObj=configDataObj,
                    xUnit="yaw",
                    errBars=errBars,
                    mkrSize=3,
                    fileName=fileName,
                    yScale="log",
                    plotTitle=(
                        "Dark Current Subtracted Data: " + scanNosString
                    ),
                )
            else:
                print(
                    "!! Please provide configDataObj if you would like to view the plot."
                )
        return

    def dataMasker(self, **kwargs):
        """
        Mask high and low gain data to remove unwanted points.

        Apply only after doing a time normalisation and dark current
        subtraction!!!

        Masking is applied to the weights, used for merging and binning later
        during processing. This allows the use of the full scan range for
        both high and low gain scans, and mitigates against USAXS scans where
        very intense peaks are seen throughout the range (i.e. gratings, etc.).

        Masking is currently accomplished by determining saturation state of 
        diodes for both high gain and low gain passes (i.e. whether you've 
        saturated the diode). If you've saturated the diode, that saturated data 
        will be removed, if not, the peak is left in. Noise level from the dark 
        current scan is used (i.e. standard deviation of the DC data) to 
        determine lower threshold limits. Currently a lower threshold is only 
        applied to high gain scans of the rear diode. You may change the default
        values for where the lower threshold is (i.e. the multiplier applied to
        the standard deviation of the dark current intensity) by passing in
        noiseMultiplier = [
            {"loGain": 5, "hiGain": False},    <- Front diode
            {"loGain": 6, "hiGain": 4}         <- Rear diode
        ]
        Integers act as multipliers for the noise threshold value, and False
        means this value is not applied (i.e. for the default the 'noisy' floor of
        the front diode for a high gain scan is kept).

        If configDataObj is passed in and plotMe is "True", data is plotted
        and plot is saved to the "/processing/Figures/" subdirectory in the
        specified dataPath.
        """
        configDataObj = kwargs.get("configDataObj", None)
        noiseMultiplier = kwargs.get(
            "noiseMultiplier",
            [{"loGain": 5, "hiGain": False}, {"loGain": 6, "hiGain": 4}],
        )
        plotMe = kwargs.get("plotMe", False)
        errBars = kwargs.get("errBars", False)

        if configDataObj is not None:
            feedback = configDataObj.feedback
        else:
            feedback = False

        print("Masking data...")
        for runNo in list(self.usaxsRuns.keys()):
            for inc in range(self.IN):
                if all(
                    (
                        self.usaxsRuns[runNo].loGain.complete[inc],
                        self.usaxsRuns[runNo].loGain.beam[inc],
                    )
                ):
                    saturated = self.usaxsRuns[runNo].loGain.saturated[inc]
                    if saturated:
                        if feedback:
                            print(
                                "Low gain scan saturated the diode, cutting out unwanted data..."
                            )
                        loGainThreshold = (
                            max(self.usaxsRuns[runNo].loGain.I[inc])
                            - (0.02 * max(self.usaxsRuns[runNo].loGain.I[inc]))
                        )  # If saturated, set threshold value slightly lower than saturation limit
                    elif not saturated:
                        if feedback:
                            print(
                                "Low gain scan did not saturate the diode, keeping high intensity data..."
                            )
                        loGainThreshold = (
                            max(self.usaxsRuns[runNo].loGain.I[inc])
                            + (0.02 * max(self.usaxsRuns[runNo].loGain.I[inc]))
                        )  # If not saturated, set threshold value slightly above than max value
                    if self.usaxsRuns[runNo].loGainDC.present:
                        loGainNoise = self.usaxsRuns[runNo].loGainDC.ISTD[inc]
                    else:
                        loGainNoise = min(self.usaxsRuns[runNo].loGain.I[inc])
                else:
                    if feedback:
                        print(
                            "!! Low gain scan either not present, incomplete or did not have beam, cannot be masked..."
                        )

                if all(
                    (
                        self.usaxsRuns[runNo].hiGain.complete[inc],
                        self.usaxsRuns[runNo].hiGain.beam[inc],
                    )
                ):
                    saturated = self.usaxsRuns[runNo].hiGain.saturated[inc]
                    if saturated:
                        if feedback:
                            print(
                                "High gain data saturated the diode, cutting out unwanted data..."
                            )
                        hiGainThreshold = (
                            max(self.usaxsRuns[runNo].hiGain.I[inc])
                            - (0.02 * max(self.usaxsRuns[runNo].hiGain.I[inc]))
                        )  # If saturated, set threshold value slightly lower than saturation limit
                    elif not saturated:
                        if feedback:
                            print(
                                "High gain data did not saturate the diode, keeping high intensity data..."
                            )
                        hiGainThreshold = (
                            max(self.usaxsRuns[runNo].hiGain.I[inc])
                            + (0.02 * max(self.usaxsRuns[runNo].hiGain.I[inc]))
                        )  # If not saturated, set threshold value slightly above than max value
                    if self.usaxsRuns[runNo].hiGainDC.present:
                        hiGainNoise = self.usaxsRuns[runNo].hiGainDC.ISTD[inc]
                    else:
                        hiGainNoise = min(self.usaxsRuns[runNo].hiGain.I[inc])
                else:
                    if feedback:
                        print(
                            "!! High gain scan either not present, incomplete or did not have beam, cannot be masked..."
                        )

                if self.usaxsRuns[runNo].loGain.present:
                    self.usaxsRuns[runNo].loGain.Mask.append(
                        np.full_like(
                            self.usaxsRuns[runNo].loGain.I[inc], True, dtype=bool
                        )
                    )
                    if noiseMultiplier[inc]["loGain"]:
                        self.usaxsRuns[runNo].loGain.Mask[inc][
                            self.usaxsRuns[runNo].loGain.I[inc]
                            < (noiseMultiplier[inc]["loGain"] * loGainNoise)
                        ] = False
                    self.usaxsRuns[runNo].loGain.Mask[inc][
                        self.usaxsRuns[runNo].loGain.I[inc] > loGainThreshold
                    ] = False
                    self.usaxsRuns[runNo].loGain.Iwt[inc][
                        ~self.usaxsRuns[runNo].loGain.Mask[inc]
                    ] = 0

                if self.usaxsRuns[runNo].hiGain.present:
                    self.usaxsRuns[runNo].hiGain.Mask.append(
                        np.full_like(
                            self.usaxsRuns[runNo].hiGain.I[inc], True, dtype=bool
                        )
                    )
                    if noiseMultiplier[inc]["hiGain"]:
                        self.usaxsRuns[runNo].hiGain.Mask[inc][
                            self.usaxsRuns[runNo].hiGain.I[inc]
                            < (noiseMultiplier[inc]["hiGain"] * hiGainNoise)
                        ] = False
                    self.usaxsRuns[runNo].hiGain.Mask[inc][
                        self.usaxsRuns[runNo].hiGain.I[inc] > hiGainThreshold
                    ] = False
                    self.usaxsRuns[runNo].hiGain.Iwt[inc][
                        ~self.usaxsRuns[runNo].hiGain.Mask[inc]
                    ] = 0

        self.processingSteps.append("Mask data")

        if plotMe:
            if configDataObj is not None:
                scanNosString = "i22"
                for scanNo in self.scanNos:
                    scanNosString += "-"
                    scanNosString += str(scanNo)
                fileName = (
                    configDataObj.dataPath
                    + "/processing/Figures/"
                    + scanNosString
                    + "_step"
                    + str(self.processingStep)
                    + "_Masking.png"
                )
                self.processingStep += 1
                scatteringDataPlotter(
                    scatteringDataObject=self,
                    configDataObj=configDataObj,
                    xUnit="yaw",
                    errBars=errBars,
                    fileName=fileName,
                    yScale="log",
                    plotTitle=(
                        "Masked Data: " + scanNosString
                    ),
                )
            else:
                print(
                    "!! Please provide configDataObj if you would like to view the plot."
                )
        return

    def diodeScale(self, **kwargs):
        """Scales front and rear diode data, either by fixed values or auto-scaled.

        You can choose to auto-scale data in a number of ways:
        - autoScale = "both" -> scales hi-gain and lo-gain independently for each 
            run
        - autoScale = "hiGain" -> calculates scaling values for hiGain, and applies
            to both hiGain and loGain passes (independently for each run)
        - autoScale = "loGain" -> calculates scaling values for loGain, and applies
            to both hiGain and loGain passes (independently for each run)
        Currently loGain scaling truncates scans to +/- 500 points around the max 
        val, while hiGain scaling uses the entire scan. This is to counteract the 
        noisy tails often present in loGain scatter. You can change these values by
        passing in something similar to scaleWidth = [500, False], where the first
        value in the list corresponds to the low gain scan, and the second the high
        gain scan. A false value will take all the points, while an integer will 
        take +/- n points around the max intensity point.

        Instead you can use default values (by passing autoScale = False, or 
        leaving it out). You may also give your own values, by passing in
        diodeScaleFactors = [
            [0.0, 1.0],           <- Values to apply to the first diode
            [0.0, 5.33]         <- Values to apply to the second diode
        ]
        Values to use take the form of a simple linear scaling:
            newIntensity = values[1] * intensity + values[0]
        """
        diodeScaleFactors = kwargs.get("diodeScaleFactors", None)
        configDataObj = kwargs.get("configDataObj", None)
        autoScale = kwargs.get("autoScale", False)
        wtScale = kwargs.get("wtScale", "log")
        scaleWidth = kwargs.get("scaleWidth",[500, False])
        plotMe = kwargs.get("plotMe", False)
        errBars = kwargs.get("errBars", False)

        if configDataObj is not None:
            feedback = configDataObj.feedback
        else:
            feedback = False

        print("Scaling diode intensitites.")

        if autoScale in ("loGain", "hiGain", "both"):
            for runNo in list(self.usaxsRuns.keys()):
                for inc in range(self.IN - 1):
                    if autoScale in ("loGain", "both"):
                        if all(
                            (
                                self.usaxsRuns[runNo].loGain.complete[inc],
                                self.usaxsRuns[runNo].loGain.beam[inc],
                            )
                        ):
                            print("Auto scaling diode intensities for loGain")
                            maxIndex = np.argmax(self.usaxsRuns[runNo].loGain.I[inc])
                            if scaleWidth[0]:
                                lowerBound = maxIndex - scaleWidth[0]
                                if lowerBound < 0:
                                    lowerBound = 0
                                upperBound = maxIndex + scaleWidth[0]
                                if upperBound > len(self.usaxsRuns[runNo].loGain.I[inc]):
                                    upperBound = len(self.usaxsRuns[runNo].loGain.I[inc])
                            elif not scaleWidth[0]:
                                lowerBound = 0
                                upperBound = len(self.usaxsRuns[runNo].loGain.I[inc])

                            ds1y = self.usaxsRuns[runNo].loGain.I[inc][
                                lowerBound:upperBound
                            ]
                            ds1e = self.usaxsRuns[runNo].loGain.Iwt[inc][
                                lowerBound:upperBound
                            ]
                            ds2y = self.usaxsRuns[runNo].loGain.I[inc + 1][
                                lowerBound:upperBound
                            ]
                            ds2e = self.usaxsRuns[runNo].loGain.Iwt[inc + 1][
                                lowerBound:upperBound
                            ]
                            divRes = ds1y / ds2y

                            sc = np.zeros(2)
                            sc[0] = 0
                            sc[1] = np.nanmean(divRes[divRes > 0])
                            if np.isnan(sc[1]):
                                sc[1] = 1

                            def csqr(sc):
                                return (ds1y - (sc[1] * ds2y) + sc[0]) / (
                                    np.sqrt((1/ds1e)**2 + (1/ds2e)**2)
                                )

                            def csqrLog(sc):
                                return np.log10(
                                    (ds1y - (sc[1] * ds2y) + sc[0])
                                    / (np.sqrt((1/ds1e)**2 + (1/ds2e)**2))
                                )

                            if wtScale == "linear":
                                sc, _ = leastsq(csqr, sc)
                            elif wtScale == "log":
                                sc, _ = leastsq(csqrLog, sc)
                            else:
                                print("!! Weighting scale not recognised...")

                            self.usaxsRuns[runNo].loGain.I[inc + 1] = (
                                sc[1] * self.usaxsRuns[runNo].loGain.I[inc + 1] + sc[0]
                            )
                            self.usaxsRuns[runNo].loGain.ISTD[inc + 1] *= abs(sc[1])
                            print(
                                f"For run number {runNo}, diode loGain scan, scaling was found to be {sc[1]} * Intensity + {sc[0]}"
                            )

                            if autoScale == "loGain":
                                if all(
                                    (
                                        self.usaxsRuns[runNo].hiGain.complete[inc],
                                        self.usaxsRuns[runNo].hiGain.beam[inc],
                                    )
                                ):
                                    print(
                                        "Auto scaling diode intensities for hiGain to values found for loGain"
                                    )
                                    self.usaxsRuns[runNo].hiGain.I[inc + 1] = (
                                        sc[1] * self.usaxsRuns[runNo].hiGain.I[inc + 1]
                                        + sc[0]
                                    )
                                    self.usaxsRuns[runNo].hiGain.ISTD[inc + 1] *= abs(
                                        sc[1]
                                    )
                                else:
                                    if feedback:
                                        print("!! hiGain scan compromised - not scaled.")
                        else:
                            if feedback:
                                print("!! loGain scan compromised, not scaled.")

                    if autoScale in ("hiGain", "both"):
                        if all(
                            (
                                self.usaxsRuns[runNo].hiGain.complete[inc],
                                self.usaxsRuns[runNo].hiGain.beam[inc],
                            )
                        ):
                            print("Auto scaling diode intensities for hiGain")
                            if scaleWidth[1]:
                                maxIndex = np.argmax(self.usaxsRuns[runNo].hiGain.I[inc])
                                lowerBound = maxIndex - scaleWidth[1]
                                if lowerBound < 0:
                                    lowerBound = 0
                                upperBound = maxIndex + scaleWidth[1]
                                if upperBound > len(self.usaxsRuns[runNo].hiGain.I[inc]):
                                    upperBound = len(self.usaxsRuns[runNo].hiGain.I[inc])
                            elif not scaleWidth[1]:
                                lowerBound = 0
                                upperBound = len(self.usaxsRuns[runNo].hiGain.I[inc])

                            ds1y = self.usaxsRuns[runNo].hiGain.I[inc][
                                lowerBound:upperBound
                            ]
                            ds1e = self.usaxsRuns[runNo].hiGain.Iwt[inc][
                                lowerBound:upperBound
                            ]
                            ds2y = self.usaxsRuns[runNo].hiGain.I[inc + 1][
                                lowerBound:upperBound
                            ]
                            ds2e = self.usaxsRuns[runNo].hiGain.Iwt[inc + 1][
                                lowerBound:upperBound
                            ]
                            divRes = ds1y / ds2y

                            sc = np.zeros(2)
                            sc[0] = 0
                            sc[1] = np.nanmean(divRes[divRes > 0])
                            if np.isnan(sc[1]):
                                sc[1] = 1

                            def csqr(sc):
                                return (ds1y - (sc[1] * ds2y) + sc[0]) / (
                                    np.sqrt((1/ds1e)**2 + (1/ds2e)**2)
                                )

                            def csqrLog(sc):
                                return np.log10(
                                    (ds1y - (sc[1] * ds2y) + sc[0])
                                    / (np.sqrt((1/ds1e)**2 + (1/ds2e)**2))
                                )

                            if wtScale == "linear":
                                sc, _ = leastsq(csqr, sc)
                            elif wtScale == "log":
                                sc, _ = leastsq(csqrLog, sc)
                            else:
                                print("!! Weighting scale not recognised...")

                            self.usaxsRuns[runNo].hiGain.I[inc + 1] = (
                                sc[1] * self.usaxsRuns[runNo].hiGain.I[inc + 1] + sc[0]
                            )
                            self.usaxsRuns[runNo].hiGain.ISTD[inc + 1] *= abs(sc[1])
                            print(
                                f"For run number {runNo}, diode hiGain scan, scaling was found to be {sc[1]} * Intensity + {sc[0]}"
                            )

                            if autoScale == "hiGain":
                                if all(
                                    (
                                        self.usaxsRuns[runNo].loGain.complete[inc],
                                        self.usaxsRuns[runNo].loGain.beam[inc],
                                    )
                                ):
                                    self.usaxsRuns[runNo].loGain.I[inc + 1] = (
                                        sc[1] * self.usaxsRuns[runNo].loGain.I[inc + 1]
                                        + sc[0]
                                    )
                                    self.usaxsRuns[runNo].loGain.ISTD[inc + 1] *= abs(
                                        sc[1]
                                    )
                                else:
                                    if feedback:
                                        print("!! loGain scan compromised - not scaled.")
                        else:
                            if feedback:
                                print("!! hiGain scan compromised, not scaled.")

        elif not autoScale:
            print("Scaling diode intensities based on fixed values.")
            if diodeScaleFactors:
                if feedback:
                    print(
                        "Scaling diode intensities based on fixed values you've specified."
                    )
                    # TODO: A check here to make sure they're the right shape...
            else:
                diodeScaleFactors = configDataObj.diodeScaleFactors
                if feedback:
                    if feedback:
                        print(
                            "Scaling diode intensities based on fixed values from configuration."
                        )

            for runNo in list(self.usaxsRuns.keys()):
                for inc in range(self.IN):
                    sc = diodeScaleFactors[inc]

                    if all(
                        (
                            self.usaxsRuns[runNo].loGain.complete[inc],
                            self.usaxsRuns[runNo].loGain.beam[inc],
                        )
                    ):
                        self.usaxsRuns[runNo].loGain.I[inc] = (
                            sc[1] * self.usaxsRuns[runNo].loGain.I[inc] + sc[0]
                        )
                        self.usaxsRuns[runNo].loGain.ISTD[inc] *= abs(sc[1])
                    else:
                        if feedback:
                            print("!! loGain scan compromised - not scaled.")

                    if all(
                        (
                            self.usaxsRuns[runNo].hiGain.complete[inc],
                            self.usaxsRuns[runNo].hiGain.beam[inc],
                        )
                    ):
                        self.usaxsRuns[runNo].hiGain.ISTD[inc] *= abs(sc[1])
                        self.usaxsRuns[runNo].hiGain.I[inc] = (
                            sc[1] * self.usaxsRuns[runNo].hiGain.I[inc] + sc[0]
                        )
                    else:
                        if feedback:
                            print("!! hiGain scan compromised - not scaled.")

        self.processingSteps.append("Match scaling")

        if plotMe:
            if configDataObj is not None:
                scanNosString = "i22"
                for scanNo in self.scanNos:
                    scanNosString += "-"
                    scanNosString += str(scanNo)
                fileName = (
                    configDataObj.dataPath
                    + "/processing/Figures/"
                    + scanNosString
                    + "_step"
                    + str(self.processingStep)
                    + "_diodeScale.png"
                )
                self.processingStep += 1
                scatteringDataPlotter(
                    scatteringDataObject=self,
                    configDataObj=configDataObj,
                    xUnit="yaw",
                    errBars=errBars,
                    fileName=fileName,
                    yScale="log",
                    plotTitle=(
                        "Diode Scaled Data: " + scanNosString
                    ),
                )
            else:
                print(
                    "!! Please provide configDataObj if you would like to view the plot."
                )
        return

    def I0Norm(self, **kwargs):
        """
        Normalise data to I0.

        Normalises intensity data to incident intensity (intensity at qbpm_2).
        If I0 is not available, or is not the correct length, nothing will be
        done. Normalise after masking!!!

        If configDataObj is passed in and plotMe is "True", data is plotted
        and plot is saved to the "/processing/Figures/" subdirectory in the
        specified dataPath.
        """
        configDataObj = kwargs.get("configDataObj", None)
        plotMe = kwargs.get("plotMe", False)
        errBars = kwargs.get("errBars", False)

        if configDataObj is not None:
            feedback = configDataObj.feedback
        else:
            feedback = False

        print("Normalising to I0...")
        for runNo in list(self.usaxsRuns.keys()):
            if all(
                (
                    self.usaxsRuns[runNo].loGain.complete
                    + self.usaxsRuns[runNo].loGain.beam
                )
            ):
                self.usaxsRuns[runNo].loGain.I0Norm()
            else:
                if feedback:
                    print(
                        "!! Low gain scan either incomplete or did not have beam, cannot be normalised to I0..."
                    )
            if all(
                (
                    self.usaxsRuns[runNo].hiGain.complete
                    + self.usaxsRuns[runNo].hiGain.beam
                )
            ):
                self.usaxsRuns[runNo].hiGain.I0Norm()
            else:
                if feedback:
                    print(
                        "!! High gain scan either incomplete or did not have beam, cannot be normalised to I0..."
                    )

        self.processingSteps.append("Normalise to I0")

        if plotMe:
            if configDataObj is not None:
                scanNosString = "i22"
                for scanNo in self.scanNos:
                    scanNosString += "-"
                    scanNosString += str(scanNo)
                fileName = (
                    configDataObj.dataPath
                    + "/processing/Figures/"
                    + scanNosString
                    + "_step"
                    + str(self.processingStep)
                    + "_I0Norm.png"
                )
                self.processingStep += 1
                scatteringDataPlotter(
                    scatteringDataObject=self,
                    configDataObj=configDataObj,
                    xUnit="yaw",
                    errBars=errBars,
                    fileName=fileName,
                    yScale="log",
                    plotTitle=(
                        "I0 Normalised Data: " + scanNosString
                    ),
                )
            else:
                print(
                    "!! Please provide configDataObj if you would like to view the plot."
                )
        return

    def dataCentrer(self, configDataObj, **kwargs):
        """
        Centre high and low gain scans around the peak of the rocking curve, IT
        (peak height of fitted gaussian) is saved.

        Centring currently accomplished by determining peak centre of low gain
        scan, and the same offset is applied to the high gain scan. To determine
        peak centre, a gaussian is fit to the low gain scan, +/- 50 data points
        around the maximum value (default). If you are forcing the centre to
        be found, it may be best to increase the widthVal (number of data
        points outside the maxVal which are used during fitting).
        If configDataObj is passed in and plotMe is "True", data is plotted
        and plot is saved to the "/processing/Figures/" subdirectory in the
        specified dataPath.
        """
        widthVal = kwargs.get("widthVal", 50)
        forceCentre = kwargs.get("forceCentre", False)
        plotMe = kwargs.get("plotMe", False)
        errBars = kwargs.get("errBars", False)
        checkFolding = kwargs.get("checkFolding", False)

        if configDataObj is not None:
            feedback = configDataObj.feedback
        else:
            feedback = False

        def centerPeak(yawVals, intensityVals, widthVal):
            try:
                maxIndex = np.argmax(intensityVals)
                scanLength = len(yawVals)

                if (maxIndex - widthVal) < 0:
                    lowerIndex = 0
                else:
                    lowerIndex = maxIndex - widthVal
                if (maxIndex + widthVal) > scanLength:
                    upperIndex = scanLength
                else:
                    upperIndex = maxIndex + widthVal

                x = yawVals[lowerIndex:upperIndex]
                y = intensityVals[lowerIndex:upperIndex]
                mean = sum(x * y) / sum(y)
                sigma = np.sqrt(sum(y * (x - mean) ** 2) / sum(y))
                # def gauss_function(b, a, b0, sigma): # Moved to usaxsToolbox.py
                #     return a*np.exp(-(b-b0)**2/(2*sigma**2))
                popt, pcov = curve_fit(gaussFunction, x, y, p0=[1, mean, sigma])
                peakCentre = popt[1]
                IT = popt[0]
                return peakCentre, IT
            except RuntimeError:
                maxIndex = np.argmax(intensityVals)
                scanLength = len(yawVals)
                widthVal = widthVal * 4

                if (maxIndex - widthVal) < 0:
                    lowerIndex = 0
                else:
                    lowerIndex = maxIndex - widthVal
                if (maxIndex + widthVal) > scanLength:
                    upperIndex = scanLength
                else:
                    upperIndex = maxIndex + widthVal

                x = yawVals[lowerIndex:upperIndex]
                y = intensityVals[lowerIndex:upperIndex]
                mean = sum(x * y) / sum(y)
                sigma = np.sqrt(sum(y * (x - mean) ** 2) / sum(y))
                # def gauss_function(b, a, b0, sigma): # Moved to usaxsToolbox.py
                #     return a*np.exp(-(b-b0)**2/(2*sigma**2))
                popt, pcov = curve_fit(gaussFunction, x, y, p0=[1, mean, sigma])
                peakCentre = popt[1]
                IT = popt[0]
                return peakCentre, IT

        print("Centering peaks...")
        for runNo in list(self.usaxsRuns.keys()):
            if all(
                (
                    self.usaxsRuns[runNo].loGain.complete
                    + self.usaxsRuns[runNo].loGain.beam
                )
            ):
                if not all((self.usaxsRuns[runNo].loGain.saturated)):
                    finiteValMask = np.isfinite(self.usaxsRuns[runNo].loGain.I[-1])
                    yaw = self.usaxsRuns[runNo].loGain.yaw[finiteValMask]
                    intensity = self.usaxsRuns[runNo].loGain.I[-1][finiteValMask]
                    peakCentre, IT = centerPeak(yaw, intensity, widthVal)
                    self.usaxsRuns[runNo].loGain.yaw = (
                        self.usaxsRuns[runNo].loGain.yaw - peakCentre
                    )
                    self.usaxsRuns[runNo].hiGain.yaw = (
                        self.usaxsRuns[runNo].hiGain.yaw - peakCentre
                    )
                    self.usaxsRuns[runNo].IT = IT
                    self.usaxsRuns[runNo].peakCentre = peakCentre
                else:
                    if feedback:
                        print(
                            "!! Low gain scan saturated the diode, cannot be centered..."
                        )
                    if forceCentre:
                        if feedback:
                            print(
                                "?? You've chosen to ignore diode saturation - centering results may be incorrect"
                            )
                        finiteValMask = np.isfinite(self.usaxsRuns[runNo].loGain.I[-1])
                        yaw = self.usaxsRuns[runNo].loGain.yaw[finiteValMask]
                        intensity = self.usaxsRuns[runNo].loGain.I[-1][finiteValMask]
                        peakCentre, IT = centerPeak(yaw, intensity, widthVal)
                        self.usaxsRuns[runNo].loGain.yaw = (
                            self.usaxsRuns[runNo].loGain.yaw - peakCentre
                        )
                        self.usaxsRuns[runNo].hiGain.yaw = (
                            self.usaxsRuns[runNo].hiGain.yaw - peakCentre
                        )
                        self.usaxsRuns[runNo].IT = IT
                        self.usaxsRuns[runNo].peakCentre = peakCentre
            else:
                if feedback:
                    print(
                        "!! Low gain scan either incomplete or did not have beam, attempting to center with high gain scan..."
                    )
                if all(
                    (
                        self.usaxsRuns[runNo].hiGain.complete
                        + self.usaxsRuns[runNo].hiGain.beam
                    )
                ):
                    if not all((self.usaxsRuns[runNo].hiGain.saturated)):
                        finiteValMask = np.isfinite(self.usaxsRuns[runNo].hiGain.I[-1])
                        yaw = self.usaxsRuns[runNo].hiGain.yaw[finiteValMask]
                        intensity = self.usaxsRuns[runNo].hiGain.I[-1][finiteValMask]
                        peakCentre, IT = centerPeak(yaw, intensity, widthVal)
                        self.usaxsRuns[runNo].hiGain.yaw = (
                            self.usaxsRuns[runNo].hiGain.yaw - peakCentre
                        )
                        self.usaxsRuns[runNo].IT = IT
                        self.usaxsRuns[runNo].peakCentre = peakCentre
                    else:
                        if feedback:
                            print(
                                "!! High gain scan saturated the diode, cannot be centered..."
                            )
                        if forceCentre:
                            if feedback:
                                print(
                                    "?? You've chosen to ignore diode saturation - centering results may be incorrect"
                                )
                            finiteValMask = np.isfinite(
                                self.usaxsRuns[runNo].hiGain.I[-1]
                            )
                            yaw = self.usaxsRuns[runNo].hiGain.yaw[finiteValMask]
                            intensity = self.usaxsRuns[runNo].hiGain.I[-1][
                                finiteValMask
                            ]
                            peakCentre, IT = centerPeak(yaw, intensity, widthVal)
                            self.usaxsRuns[runNo].hiGain.yaw = (
                                self.usaxsRuns[runNo].hiGain.yaw - peakCentre
                            )
                            self.usaxsRuns[runNo].IT = IT
                            self.usaxsRuns[runNo].peakCentre = peakCentre

                else:
                    if feedback:
                        print(
                            "!! High gain scan either incomplete or did not have beam, cannot be centered..."
                        )

        self.processingSteps.append("Center data")

        if checkFolding:
            if configDataObj is not None:
                scanNosString = "i22"
                for scanNo in self.scanNos:
                    scanNosString += "-"
                    scanNosString += str(scanNo)
                fileName = (
                    configDataObj.dataPath
                    + "/processing/Figures/"
                    + scanNosString
                    + "_step"
                    + str(self.processingStep)
                    + "_FoldingCheck.png"
                )
                self.processingStep += 1
                foldedDataPlotter(
                    scatteringDataObject=self,
                    configDataObj=configDataObj,
                    xUnit="yaw",
                    fileName=fileName,
                    # xScale="log",
                    yScale="log",
                    plotTitle=(
                        "Centred Data: " + scanNosString
                    ),
                )
            else:
                print(
                    "!! Please provide configDataObj if you would like to view the plot."
                )

        if plotMe:
            if configDataObj is not None:
                scanNosString = "i22"
                for scanNo in self.scanNos:
                    scanNosString += "-"
                    scanNosString += str(scanNo)
                fileName = (
                    configDataObj.dataPath
                    + "/processing/Figures/"
                    + scanNosString
                    + "_step"
                    + str(self.processingStep)
                    + "_Centering.png"
                )
                self.processingStep += 1
                scatteringDataPlotter(
                    scatteringDataObject=self,
                    configDataObj=configDataObj,
                    xUnit="yaw",
                    errBars=errBars,
                    fileName=fileName,
                    yScale="log",
                    plotTitle=(
                        "Centred Data: " + scanNosString
                    ),
                )
            else:
                print(
                    "!! Please provide configDataObj if you would like to view the plot."
                )
        return

    def transmissionNorm(self, **kwargs):
        """
        Normalise data to transmission value (i.e. peak intensity of rocking
        curve).

        Normalises intensity (I) data to transmission value found in centered
        data.
        If configDataObj is passed in and plotMe is "True", data is plotted
        and plot is saved to the "/processing/Figures/" subdirectory in the
        specified dataPath.
        """
        configDataObj = kwargs.get("configDataObj", None)
        plotMe = kwargs.get("plotMe", False)
        errBars = kwargs.get("errBars", False)

        if configDataObj is not None:
            feedback = configDataObj.feedback
        else:
            feedback = False

        print("Normalising to IT...")
        for runNo in list(self.usaxsRuns.keys()):
            try:
                self.usaxsRuns[runNo].IT.size
                IT = self.usaxsRuns[runNo].IT
            except:
                IT = 1
                if feedback:
                    print("!! No IT value was found, using 1...")
            if all(
                (
                    self.usaxsRuns[runNo].loGain.complete
                    + self.usaxsRuns[runNo].loGain.beam
                )
            ):
                self.usaxsRuns[runNo].loGain.I /= IT
                self.usaxsRuns[runNo].loGain.ISTD /= IT
            if all(
                (
                    self.usaxsRuns[runNo].hiGain.complete
                    + self.usaxsRuns[runNo].hiGain.beam
                )
            ):
                self.usaxsRuns[runNo].hiGain.I /= IT
                self.usaxsRuns[runNo].hiGain.ISTD /= IT

        self.processingSteps.append("Normalise to transmission")

        if plotMe:
            if configDataObj is not None:
                scanNosString = "i22"
                for scanNo in self.scanNos:
                    scanNosString += "-"
                    scanNosString += str(scanNo)
                fileName = (
                    configDataObj.dataPath
                    + "/processing/Figures/"
                    + scanNosString
                    + "_step"
                    + str(self.processingStep)
                    + "_ITNorm.png"
                )
                self.processingStep += 1
                scatteringDataPlotter(
                    scatteringDataObject=self,
                    configDataObj=configDataObj,
                    xUnit="yaw",
                    errBars=errBars,
                    fileName=fileName,
                    yScale="log",
                    plotTitle=(
                        "Transmission Normalised Data: " + scanNosString
                    ),
                )
            else:
                print(
                    "!! Please provide configDataObj if you would like to view the plot."
                )
        return

    def dataSorter(self, **kwargs):
        """
        'Fold' and sort data based on yaw values for high and low gain scans.

        You have the option to choose which side(s) is(are) retained for
        further processing: both, positive, or negative. If oxidation is noted,
        it is recommended you retain the positive side. To preview folding
        results by side, please use argument 'checkFolding=True' when
        centering data.

        If configDataObj is passed in and plotMe is "True", data is plotted
        and plot is saved to the "/processing/Figures/" subdirectory in the
        specified dataPath.
        """
        side = kwargs.get("side", "both")
        configDataObj = kwargs.get("configDataObj", None)
        plotMe = kwargs.get("plotMe", False)
        errBars = kwargs.get("errBars", False)

        if configDataObj is not None:
            feedback = configDataObj.feedback
        else:
            feedback = False

        print("Folding Data...")
        for runNo in list(self.usaxsRuns.keys()):
            if all(
                (
                    self.usaxsRuns[runNo].loGain.complete
                    + self.usaxsRuns[runNo].loGain.beam
                )
            ):
                self.usaxsRuns[runNo].loGain.dataSorter(side=side)
                if feedback:
                    print("Sorted low gain data")
            if all(
                (
                    self.usaxsRuns[runNo].hiGain.complete
                    + self.usaxsRuns[runNo].hiGain.beam
                )
            ):
                self.usaxsRuns[runNo].hiGain.dataSorter(side=side)
                if feedback:
                    print("Sorted high gain data")

        self.processingSteps.append("Fold data - " + side)

        if plotMe:
            if configDataObj is not None:
                scanNosString = "i22"
                for scanNo in self.scanNos:
                    scanNosString += "-"
                    scanNosString += str(scanNo)
                fileName = (
                    configDataObj.dataPath
                    + "/processing/Figures/"
                    + scanNosString
                    + "_step"
                    + str(self.processingStep)
                    + "_Folded.png"
                )
                self.processingStep += 1
                scatteringDataPlotter(
                    scatteringDataObject=self,
                    configDataObj=configDataObj,
                    xUnit="yaw",
                    errBars=errBars,
                    fileName=fileName,
                    # xScale="log",
                    yScale="log",
                    plotTitle=(
                        "Folded Data: " + scanNosString
                    ),
                )
            else:
                print(
                    "!! Please provide configDataObj if you would like to view the plot."
                )
        return

    def qConverter(self, configDataObj, **kwargs):
        """
        Convert yaw values to q for high and low gain scans. Units for q are
        specified in the configDataObj.

        If plotMe is "True", data is plotted and plot is saved to the
        "/processing/Figures/" subdirectory in the specified dataPath.
        """
        plotMe = kwargs.get("plotMe", False)
        errBars = kwargs.get("errBars", False)

        if configDataObj is not None:
            feedback = configDataObj.feedback
        else:
            feedback = False

        print("Converting from yaw to q...")
        for runNo in list(self.usaxsRuns.keys()):
            if all(
                (
                    self.usaxsRuns[runNo].loGain.complete
                    + self.usaxsRuns[runNo].loGain.beam
                )
            ):
                self.usaxsRuns[runNo].loGain.qConverter(configDataObj)
                if feedback:
                    print("Converted low gain yaw to q")
            if all(
                (
                    self.usaxsRuns[runNo].hiGain.complete
                    + self.usaxsRuns[runNo].hiGain.beam
                )
            ):
                self.usaxsRuns[runNo].hiGain.qConverter(configDataObj)
                if feedback:
                    print("Converted high gain yaw to q")
        self.processingSteps.append("Convert to q")

        if plotMe:
            if configDataObj is not None:
                scanNosString = "i22"
                for scanNo in self.scanNos:
                    scanNosString += "-"
                    scanNosString += str(scanNo)
                fileName = (
                    configDataObj.dataPath
                    + "/processing/Figures/"
                    + scanNosString
                    + "_step"
                    + str(self.processingStep)
                    + "_qConversion.png"
                )
                self.processingStep += 1
                scatteringDataPlotter(
                    scatteringDataObject=self,
                    configDataObj=configDataObj,
                    xUnit="Q",
                    errBars=errBars,
                    fileName=fileName,
                    xScale="log",
                    yScale="log",
                    plotTitle=(
                        "q Converted Data: " + scanNosString
                    ),
                )
            else:
                print(
                    "!! Please provide configDataObj if you would like to view the plot."
                )
        return

    def dataMergeAndBin(self, configDataObj, **kwargs):
        """
        Merge and bin all USAXS scans for one sample.

        All data is merged and binned, with statistics (i.e. ISTD, QSTD, etc.)
        weighted based on statistics from tetramm and masking.

        If configDataObj is passed in and plotMe is "True", data is plotted
        and plot is saved to the "/processing/Figures/" subdirectory in the
        specified dataPath. Determined scale factor can be saved to
        configDataObj for future use.
        """
        plotMe = kwargs.get("plotMe", False)
        errBars = kwargs.get("errBars", True)
        xUnit = kwargs.get("xUnit", "Q")

        if configDataObj is not None:
            feedback = configDataObj.feedback
        else:
            feedback = False

        print("Merging and binning all usaxs runs...")

        allX = np.array([])
        allI = np.array([])
        allIwt = np.array([])
        allISTD = np.array([])
        scanNosUsed = list()

        for runNo in list(self.usaxsRuns.keys()):
            if all(
                (
                    self.usaxsRuns[runNo].loGain.complete + self.usaxsRuns[runNo].loGain.beam
                )
            ):
                for inc in range(self.IN):
                    if xUnit == "Q":
                        allX = np.concatenate((allX, self.usaxsRuns[runNo].loGain.Q))
                    elif xUnit == "yaw":
                        allX = np.concatenate((allX, self.usaxsRuns[runNo].loGain.yaw))
                    allI = np.concatenate((allI, self.usaxsRuns[runNo].loGain.I[inc]))
                    allIwt = np.concatenate(
                        (allIwt, self.usaxsRuns[runNo].loGain.Iwt[inc])
                    )
                    allISTD = np.concatenate(
                        (allISTD, self.usaxsRuns[runNo].loGain.ISTD[inc])
                    )
                scanNosUsed.append(self.usaxsRuns[runNo].loGain.scanNo)

            if all(
                (
                    self.usaxsRuns[runNo].hiGain.complete + self.usaxsRuns[runNo].hiGain.beam
                )
            ):
                for inc in range(self.IN):
                    if xUnit == "Q":
                        allX = np.concatenate((allX, self.usaxsRuns[runNo].hiGain.Q))
                    elif xUnit == "yaw":
                        allX = np.concatenate((allX, self.usaxsRuns[runNo].hiGain.yaw))
                    allI = np.concatenate((allI, self.usaxsRuns[runNo].hiGain.I[inc]))
                    allIwt = np.concatenate(
                        (allIwt, self.usaxsRuns[runNo].hiGain.Iwt[inc])
                    )
                    allISTD = np.concatenate(
                        (allISTD, self.usaxsRuns[runNo].hiGain.ISTD[inc])
                    )
                scanNosUsed.append(self.usaxsRuns[runNo].hiGain.scanNo)

        # Sort high and low gain scans
        inds = allX.argsort()
        allX = allX[inds]
        allI = allI[inds]
        allISTD = allISTD[inds]
        allIwt = allIwt[inds]

        # Init nan-ful arrays for binning
        if xUnit == "Q":
            binEdges = configDataObj.binEdges
            xArray = configDataObj.qArray
        elif xUnit == "yaw":
            binEdges = np.linspace(min(allX), max(allX), (configDataObj.nBins + 2))
            xArray = np.linspace(min(allX), max(allX), (configDataObj.nBins + 2))

        IData = np.full(len(binEdges) - 1, np.nan)
        X = np.full(len(binEdges) - 1, np.nan)
        ISTD = np.full(len(binEdges) - 1, np.nan)
        ISEM = np.full(len(binEdges) - 1, np.nan)
        XSTD = np.full(len(binEdges) - 1, np.nan)
        XSEM = np.full(len(binEdges) - 1, np.nan)
        N = np.full(len(binEdges) - 1, np.nan)

        edgeIndices = np.searchsorted(allX, binEdges)

        for binN in range(len(binEdges) - 1):
            lowerIndex, upperIndex = edgeIndices[binN], edgeIndices[binN + 1]
            rangeLen = upperIndex - lowerIndex
            if rangeLen == 0:  # Nothing to bin
                N[binN] = 0
            elif rangeLen == 1:  # One data point to bin
                N[binN] = 1
                IData[binN] = float(allI[lowerIndex:upperIndex][0])
                X[binN] = float(allX[lowerIndex:upperIndex][0])
            else:  # Multiple data points in bin
                N[binN] = np.count_nonzero(allIwt[lowerIndex:upperIndex])
                DSI = DescrStatsW(
                    allI[lowerIndex:upperIndex], weights=allIwt[lowerIndex:upperIndex]
                )
                DSX = DescrStatsW(
                    allX[lowerIndex:upperIndex], weights=allIwt[lowerIndex:upperIndex]
                )

                IData[binN] = DSI.mean
                X[binN] = DSX.mean

                ISTD[binN] = DSI.std  # TODO: Propagate original tetramm errors through
                ISEM[binN] = DSI.std * np.sqrt(
                    (allIwt[lowerIndex:upperIndex] ** 2).sum()
                    / (allIwt[lowerIndex:upperIndex].sum()) ** 2
                )

                XSTD[binN] = DSX.std
                XSEM[binN] = DSX.std * np.sqrt(
                    (allIwt[lowerIndex:upperIndex] ** 2).sum()
                    / (allIwt[lowerIndex:upperIndex].sum()) ** 2
                )

        # Remove nans
        nanMask = ~np.isnan(IData)
        IData = IData[nanMask]
        ISTD = ISTD[nanMask]
        ISEM = ISEM[nanMask]
        X = X[nanMask]
        XSTD = XSTD[nanMask]
        XSEM = XSEM[nanMask]
        N = N[nanMask]
        interpX = xArray

        # Interpolate onto the original requested qArray
        fI = interp1d(X, IData, kind="linear", bounds_error=False)
        interpI = fI(interpX)
        FISTD = interp1d(X, ISTD, kind="linear", bounds_error=False)
        interpISTD = FISTD(interpX)
        FISEM = interp1d(X, ISTD, kind="linear", bounds_error=False)
        interpISEM = FISEM(interpX)
        FXSTD = interp1d(X, XSTD, kind="linear", bounds_error=False)
        interpXSTD = FXSTD(interpX)
        FXSEM = interp1d(X, XSEM, kind="linear", bounds_error=False)
        interpXSEM = FXSEM(interpX)
        FBinnedN = interp1d(X, N, kind="linear", bounds_error=False)
        interpBinnedN = FBinnedN(interpX)

        nanMask = ~np.isnan(interpI)
        interpI = interpI[nanMask]
        interpISTD = interpISTD[nanMask]
        interpISEM = interpISEM[nanMask]
        interpX = interpX[nanMask]
        interpXSTD = interpXSTD[nanMask]
        interpXSEM = interpXSEM[nanMask]
        interpBinnedN = interpBinnedN[nanMask]

        if xUnit == "Q":
            self.mData.setQ(interpX, interpXSTD, interpXSEM)
        elif xUnit == "yaw":
            self.mData.setYaw(interpX, interpXSTD, interpXSEM)
        self.mData.setI(interpI, interpISTD, interpISEM)
        self.mData.setDataScanNos(scanNosUsed)
        self.mData.binnedN = interpBinnedN

        self.processingSteps.append("Merge and bin data")

        if plotMe:
            if configDataObj is not None:
                scanNosString = "i22"
                for scanNo in self.scanNos:
                    scanNosString += "-"
                    scanNosString += str(scanNo)
                fileName = (
                    configDataObj.dataPath
                    + "/processing/Figures/"
                    + scanNosString
                    + "_step"
                    + str(self.processingStep)
                    + "_MergingAndBinning.png"
                )
                self.processingStep += 1
                if xUnit == "Q":
                    xScale = "log"
                elif xUnit == "yaw":
                    xScale = "linear"
                mergedDataPlotter(
                    scatteringDataObject=self,
                    configDataObj=configDataObj,
                    xUnit=xUnit,
                    errBars=errBars,
                    fileName=fileName,
                    xScale=xScale,
                    yScale="log",
                    plotTitle=(
                        "Merged and Binned Data: " + scanNosString
                    ),
                )
            else:
                print(
                    "!! Please provide configDataObj if you would like to view the plot."
                )
        return

    def backgroundSubtracter(
        self, backgroundScatteringDataObj, configDataObj, **kwargs
    ):
        """
        Subtract background (merged data) from sample (merged data).

        Background scan numbers are saved to self.mData.
        If configDataObj is passed in and plotMe is "True", data is plotted
        and plot is saved to the "/processing/Figures/" subdirectory in the
        specified dataPath.
        """
        plotMe = kwargs.get("plotMe", False)
        errBars = kwargs.get("errBars", False)

        if configDataObj is not None:
            feedback = configDataObj.feedback
        else:
            feedback = False

        print("Subtracting background...")

        # TODO: Displaced volume correction?

        if (
            len(backgroundScatteringDataObj.mData.I) == 0
            or len(backgroundScatteringDataObj.mData.Q) == 0
        ):
            print(
                "!! You must merge and bin high gain and low gain datasets to enable background subtraction."
            )
        else:
            if len(backgroundScatteringDataObj.mData.Q) == len(self.mData.Q):
                INew = self.mData.I - backgroundScatteringDataObj.mData.I
                self.mData.I = INew
                self.mData.ISTD = (
                    self.mData.ISTD + backgroundScatteringDataObj.mData.ISTD
                )
                # TODO: ISEM
            else:  # NANs got removed at some point - put them back onto the normal grid
                bgI = np.full(len(configDataObj.qArray), np.nan)
                bgISTD = np.full(len(configDataObj.qArray), np.nan)
                bgQSTD = np.full(len(configDataObj.qArray), np.nan)
                bgBinnedN = np.full(len(configDataObj.qArray), np.nan)
                dataI = np.full(len(configDataObj.qArray), np.nan)
                dataISTD = np.full(len(configDataObj.qArray), np.nan)
                dataQ = configDataObj.qArray
                dataQSTD = np.full(len(configDataObj.qArray), np.nan)
                dataBinnedN = np.full(len(configDataObj.qArray), np.nan)
                for ii in range(len(backgroundScatteringDataObj.mData.Q)):
                    index = np.argwhere(
                        configDataObj.qArray == backgroundScatteringDataObj.mData.Q[ii]
                    )
                    bgI[index] = backgroundScatteringDataObj.mData.I[ii]
                    bgISTD[index] = backgroundScatteringDataObj.mData.ISTD[ii]
                    bgQSTD[index] = backgroundScatteringDataObj.mData.QSTD[ii]
                    bgBinnedN[index] = backgroundScatteringDataObj.mData.binnedN[ii]
                for ii in range(len(self.mData.Q)):
                    index = np.argwhere(configDataObj.qArray == self.mData.Q[ii])
                    dataI[index] = self.mData.I[ii]
                    dataISTD[index] = self.mData.ISTD[ii]
                    dataQSTD[index] = self.mData.QSTD[ii]
                    dataBinnedN[index] = self.mData.binnedN[ii]

                self.mData.I = dataI - bgI
                self.mData.ISTD = dataISTD + bgISTD
                self.mData.Q = dataQ
                self.mData.QSTD = dataQSTD
                self.mData.binnedN = dataBinnedN
                # TODO: ISEM

        self.mData.setBGScanNos(backgroundScatteringDataObj.mData.dataScanNos)

        self.processingSteps.append("Subtract background")

        if plotMe:
            if configDataObj is not None:
                scanNosString = "i22"
                for scanNo in self.scanNos:
                    scanNosString += "-"
                    scanNosString += str(scanNo)
                fileName = (
                    configDataObj.dataPath
                    + "/processing/Figures/"
                    + scanNosString
                    + "_step"
                    + str(self.processingStep)
                    + "_BackgroundSub.png"
                )
                self.processingStep += 1
                mergedDataPlotter(
                    scatteringDataObject=self,
                    configDataObj=configDataObj,
                    xUnit="Q",
                    errBars=errBars,
                    fileName=fileName,
                    xScale="log",
                    yScale="log",
                    plotTitle=(
                        "Background Subtracted Data: " + scanNosString
                    ),
                )
            else:
                print(
                    "!! Please provide configDataObj if you would like to view the plot."
                )
        return

    def processedDataMultiplier(self, configDataObj, **kwargs):
        """
        Multiply merged data by a scalar sepcified in configDataObj.
        """
        plotMe = kwargs.get("plotMe", False)
        errBars = kwargs.get("errBars", False)

        if configDataObj is not None:
            feedback = configDataObj.feedback
        else:
            feedback = False

        if feedback:
            print(
                "Multiplying by a set scalar ("
                + str(configDataObj.scalarMultipler)
                + ")."
            )
        self.mData.I *= configDataObj.scalarMultipler
        self.mData.ISTD *= configDataObj.scalarMultipler
        self.mData.ISEM *= configDataObj.scalarMultipler

        self.processingSteps.append("Fixed multiplier")

        return

    def qCalibrationGuess(self, configDataObj, mulVal=1, **kwargs):
        """
        For use only with SAXS and de-smeared USAXS data. Use with caution.
        """
        # TODO: Fix this if we need it...
        plotMe = kwargs.get("plotMe", True)

        if len(self.mData.desmearedQ) >= 1:
            ds1q = self.mData.saxsQ
            ds1y = self.mData.saxsI
            ds1e = self.mData.saxsIerr
            ds2q = self.mData.desmearedQ
            ds2y = self.mData.desmearedI
            ds2e = self.mData.desmearedIerr
        else:
            print("You need to import desmeared USAXS data to use this function...")

        minQ = max((min(ds1q), min(ds2q)))
        maxQ = min((max(ds1q), max(ds2q)))

        ds1mask = np.full(len(ds1q), True)
        ds1mask[ds1q < minQ] = False
        ds1mask[ds1q > maxQ] = False

        ds2mask = np.full(len(ds2q), True)
        ds2mask[ds2q < minQ] = False
        ds2mask[ds2q > maxQ] = False

        ds1q = ds1q[ds1mask]
        ds1y = ds1y[ds1mask]
        ds1e = ds1e[ds1mask]

        ds2q_original = ds2q  # [ds2mask]
        ds2y_original = ds2y  # [ds2mask]
        ds2e_original = ds2e  # [ds2mask]

        yaw_original = (
            np.arcsin(
                ds2q_original
                * configDataObj.outputConverter
                * ((sp.h * sp.c) / (4 * np.pi * configDataObj.energy))
            )
            * 2
        )
        ds2q_guessed = (
            (4 * np.pi * configDataObj.energy)
            / (sp.h * sp.c)
            * np.sin((yaw_original * mulVal) / 2)
            / configDataObj.outputConverter
        )

        qInterp = ds1q

        ds1yI = interp1d(ds1q, ds1y, kind="linear", bounds_error=False)
        ds1y = ds1yI(qInterp)
        ds1eI = interp1d(ds1q, ds1e, kind="linear", bounds_error=False)
        ds1e = ds1eI(qInterp)
        ds1q = qInterp

        def csqr(ds1y, ds1e, ds2y, ds2e):
            return (ds1y - ds2y) / (np.sqrt(ds1e**2 + ds2e**2))

        def csqrLog(ds1y, ds1e, ds2y, ds2e):
            return np.log10((ds1y - ds2y) / (np.sqrt(ds1e**2 + ds2e**2)))

        ds2yI = interp1d(
            ds2q_guessed * mulVal, ds2y_original, kind="linear", bounds_error=False
        )
        ds2y = ds2yI(qInterp)
        ds2eI = interp1d(
            ds2q_guessed * mulVal, ds2e_original, kind="linear", bounds_error=False
        )
        ds2e = ds2eI(qInterp)
        ds2q = qInterp

        residuals = csqr(ds1y, ds1e, ds2y, ds2e)

        if plotMe:
            fig, axs = plt.subplots(nrows=2, sharex=True)
            axs[0].errorbar(
                ds2q_original,
                ds2y_original,
                yerr=ds2e_original,
                ms=2,
                fmt="o",
                label="Original USAXS",
                alpha=0.1,
            )
            axs[0].errorbar(
                ds2q_guessed,
                ds2y_original,
                yerr=ds2e_original,
                ms=2,
                fmt="o",
                label="Modified USAXS",
                alpha=0.25,
            )
            axs[0].errorbar(
                self.mData.saxsQ,
                self.mData.saxsI,
                yerr=self.mData.saxsIerr,
                ms=5,
                fmt="o",
                label="SAXS",
                alpha=0.75,
            )
            axs[0].set_title("Multiplier Value: " + str(mulVal))
            axs[0].set_yscale("log")
            axs[0].set_ylim(min(ds1y) * 0.5, max(ds1y) * 1.5)
            axs[0].legend()

            axs[1].scatter(
                qInterp,
                residuals,
                label=(
                    "mean residual = "
                    + str(np.round(np.nanmean(abs(residuals)), decimals=4))
                ),
                alpha=0.25,
            )
            axs[1].set_xlim(minQ, maxQ)
            axs[1].legend()

            plt.show()
            plt.savefig(
                (
                    configDataObj.dataPath
                    + "/processing/usaxsQCalibration/"
                    + str(np.round(mulVal, 4))
                    + ".png"
                )
            )
            plt.close()
        return residuals
