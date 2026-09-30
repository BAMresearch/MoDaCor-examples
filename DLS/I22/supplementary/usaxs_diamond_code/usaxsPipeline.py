"""
To run this pipeline, you will need the supporting files in place:
- usaxsDataClasses
- usaxsIO
- usaxsPlotter
- usaxsToolbox

In the correct virtual environment, type "python -W ignore <pathToPipeline>/usaxsPipeline-generic.py"
"""

from usaxsDataClasses import configDataObj
from usaxsIO import (
    divvyScanNos,
    backgroundDataGrabber,
    sampleDataGrabber,
    xyeWriter,
    nexusWriter,
)
import numpy as np

###############################################################################
# Set up data objects                                                         #
###############################################################################

subtractBackground = True
lazyLoad = True
saveBackgroundFiles = ["nexus"] # Options: "nexus" & "xye"
saveSampleFiles = ["nexus"]     # Options: "nexus" & "xye"

backgroundStartNo = 977723 # Empty furnace
backgroundEndNo = 977727

sampleStartNo = 978548 # 966603
sampleEndNo = 978561 # 966612

# Shared processing setup
treatFirstYawValue = True # Manually set first yaw value to np.nan - if acquisition is a bit janky
maskingNoiseMultiplier = [{"loGain": 5, "hiGain": False}, {"loGain": 7, "hiGain": 7}]
diodeAutoScale = "both" # Scan to use for diode scaling: "loGain", "hiGain", "both", False
scaleWidth = [100, 100] # Range(s) to use for diode scaling [loGain, hiGain]
forceCentre = True # Force centering if diodes are determined to be saturated: True, False
side = "both" # 'Side' of yaw data to use for processing: "positive", "negative", "both"
scalarMultiplier = True # Multipy final data for visualisation purposes: True, False
plotMe = False # Save all plots from processing pipeline: True, False

###############################################################################
# Set up configuration object                                                 #
###############################################################################

config = configDataObj()
config.setDataPath("/dls/i22/data/2026/sm43108-1")

config.setEnergy(14)
config.setFWHM(11.33)
config.calcQMin()

config.setRequestedFrames((2000.0, 20000.0))  # Dark current then USAXS scans
config.setDarkCurrentOffset = [[3.02e-9,2e-11,],[8e-10,8e-10,]]
config.setNBins(1024)
config.calcQArray()
config.calcBinEdges()

config.setQCalibrationMultipler(0.9999935807521158)
config.setScalarMultipler(5e8)

config.feedback = False  # Make processing verbose / not verbose

###############################################################################
# Background and Sample Processing Order should be:   -\                      #
# -> setScanDetails()                                  |                      #
# -> dataShaper()                                      |                      #
# -> scanChecker()                                     |                      #
# -> timeNorm()                                        |                      #
# -> darkCurrentSubtraction()                          |                      #
# -> dataMasker()                                      |- -> Shared steps     #
# -> diodeScale()                                      |                      #
# -> I0Norm()                                          |                      #
# -> dataCenterer()                                    |                      #
# -> transmissionNorm()                                |                      #
# -> dataSorter()                                      |                      #
# -> qConverter()                                      |                      #
# -> dataMergeAndBin()                                -/                      #
# -> backgroundSubtracter()                                                   #
# -> processedDataMultiplier()                                                #
###############################################################################

###############################################################################
# Set up shared processing steps                                              #
###############################################################################

def sharedProcessing(scatteringData,config=config):
    for name in list(scatteringData.keys()):
        print(f"\nNow processing sample: {name}")
        if scatteringData[name].usaxsPresent:
            scatteringData[name].setScanDetails(configDataObj=config)
            scatteringData[name].dataShaper(configDataObj=config, plotMe=plotMe, treatFirstValue=treatFirstYawValue)
            if scatteringData[name].scanChecker(config):
                scatteringData[name].timeNorm(configDataObj=config, plotMe=plotMe)
                scatteringData[name].darkCurrentSubtraction(
                    configDataObj=config, plotMe=plotMe
                )
                scatteringData[name].dataMasker(configDataObj=config, plotMe=plotMe, noiseMultiplier=maskingNoiseMultiplier)
                scatteringData[name].diodeScale(
                    configDataObj=config, autoScale=diodeAutoScale, scaleWidth=scaleWidth, plotMe=plotMe
                )
                scatteringData[name].I0Norm(configDataObj=config, plotMe=plotMe)
                scatteringData[name].dataCentrer(
                    forceCentre=forceCentre, configDataObj=config, plotMe=plotMe, checkFolding=True
                )
                scatteringData[name].transmissionNorm(configDataObj=config, plotMe=plotMe)
                scatteringData[name].dataSorter(configDataObj=config, plotMe=plotMe, side=side)
                scatteringData[name].qConverter(configDataObj=config, plotMe=plotMe)
                scatteringData[name].dataMergeAndBin(
                    configDataObj=config, plotMe=True, errBars="STD"
                )
                scatteringData[name].sharedProcessingDone = True
            else:
                print(f"\nSample {name} was incomplete, could not be processed...")
        else:
             print(f"\nSample {name} held no USAXS data, could not be processed...")
    return scatteringData

###############################################################################
# Background Processing                                                       #
###############################################################################
if subtractBackground:
    if lazyLoad:
        scansNos = divvyScanNos(
            np.arange(backgroundStartNo, backgroundEndNo + 1, 1), configDataObj=config
        )  
        print(f"\nFound {len(scansNos)} background sample(s) to process.")
    elif not lazyLoad:
        scansNos = [np.arange(backgroundStartNo, backgroundEndNo + 1, 1)]

    print("\n----------------------------------------------------------------------\n")

    config.setBackgroundNos(scansNos[0])

    print("\n----------------------------------------------------------------------\n")

    bgData = {}
    bgData = backgroundDataGrabber(configDataObj=config)

    bgData = sharedProcessing(bgData)
    for name in list(bgData.keys()):
        if bgData[name].sharedProcessingDone:
            if "nexus" in saveBackgroundFiles:
                for name in list(bgData.keys()):
                    nexusWriter("usaxsData", bgData[name], config)

            if "xye" in saveBackgroundFiles:
                for name in list(bgData.keys()):
                    xyeWriter("usaxsData", bgData[name], config)
        else:
            print(f"Shared processing routine was not completed for {name}... \nFurther processing cannot be done for this sample.")

    if len(bgData) > 1:
        print("More than one USAXS collection was given as a background - only the first scan will be used for further processing.")

    print("\n----------------------------------------------------------------------\n")

###############################################################################
# Sample Processing                                                           #
###############################################################################
if lazyLoad:
    scansNos = divvyScanNos(
        np.arange(sampleStartNo, sampleEndNo + 1, 1), configDataObj=config
    )  
    print(f"\nFound {len(scansNos)} sample(s) to process.")

elif not lazyLoad:
    scansNos = [np.arange(sampleStartNo, sampleEndNo + 1, 1)]

for ii in range(len(scansNos)):
    config.setScanNos(scansNos[ii])
    plotMe = True

    print("\n----------------------------------------------------------------------\n")

    scanData = {}
    scanData = sampleDataGrabber(configDataObj=config)

    scanData = sharedProcessing(scanData)
    for name in list(scanData.keys()):
        if scanData[name].sharedProcessingDone:
            if subtractBackground:
                scanData[name].backgroundSubtracter(bgData[list(bgData.keys())[0]],configDataObj=config,plotMe=plotMe,errBars="STD")
            if scalarMultiplier:
                scanData[name].processedDataMultiplier(configDataObj=config, plotMe=plotMe, errBars="STD")

            if "nexus" in saveSampleFiles:
                nexusWriter("usaxsData-bgsbu", scanData[name], config)

            if "xye" in saveSampleFiles:
                    xyeWriter("usaxsData-bgsbu", scanData[name], config)
        else:
            print(f"Shared processing routine was not completed for {name}... \nFurther processing cannot be done for this sample.")

print("\n----------------------------------------------------------------------\n")
print("All done processing!")
