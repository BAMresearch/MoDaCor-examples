import os
import numpy as np
import scipy.constants as sp

def directoryChecker(fullPath):
    """
    Check if a directory exists, create if it doesn't.
    """
    os.makedirs(fullPath, mode=0o777, exist_ok=True)
    
def fileChecker(fullPath):
    """
    Check if a file exists, delete it if it does.
    """
    try:
        os.remove(fullPath)
    except OSError:
        pass      

def findNearest(array,value):
    """
    Find index of nearest value in an array.
    """
    array = np.asarray(array)
    idx = (np.abs(array-value)).argmin()
    return idx

def gaussFunction(b, a, b0, sigma):
    return a*np.exp(-(b-b0)**2/(2*sigma**2))

def saturationChecker(intensities, saturationWidth = 15, intensityVariation=0.02):
    """
    Check if you've saturated your diode
    """
    if np.sum([intensities > (max(intensities)*(1-intensityVariation))],) > saturationWidth: # If you have more than 15 scan points which same value +/- 2 %, it's saturated
        return True
    else:
        return False

def tetrammShaper(tetrammData, channel=None):
    """
    Extract the correct provided channel and average points for provided
    tetramm data. 
    
    If no channel is provided, assume only one channel is passed. 
    """
    if channel:
        A = []
        sA = []
        wA = []
        for dataset in tetrammData:
            if len(dataset) >= 1:
                tetrammDims = dataset.ndim
                if tetrammDims == 3:
                    tmpA = np.average(dataset[:,:,channel], axis=1)
                    tmpsA = np.std(dataset[:,:,channel], axis=1)
                    # wA = A/(sA**2)
                    tmpwA = 1/(tmpsA**2) # Cannot weight by A, as this artifically supresses points near 0 
                    A.append(tmpA)
                    sA.append(tmpsA)
                    wA.append(abs(tmpwA))
                elif tetrammDims == 1:
                    tmpA = dataset
                    tmpsA = 1/tmpA # Have a guess so there's an array there...
                    tmpwA = np.ones_like(tmpA)
                    A.append(tmpA)
                    sA.append(tmpsA)
                    wA.append(abs(tmpwA))
                else:
                    print("Unexpected shape of tetramm data.")
    else:
        A = []
        sA = []
        wA = []
        for dataset in tetrammData:
            if len(dataset) >= 1:
                tetrammDims = dataset.ndim
                if tetrammDims == 2:
                    tmpA = np.average(dataset, axis=1)
                    tmpsA = np.std(dataset, axis=1)
                    # wA = A/(sA**2)
                    tmpwA = 1/(tmpsA**2) # Cannot weight by A, as this artifically supresses points near 0 
                    A.append(tmpA)
                    sA.append(tmpsA)
                    wA.append(abs(tmpwA))
                elif tetrammDims == 1:
                    tmpA = dataset
                    tmpsA = 1/tmpA # Have a guess so there's an array there...
                    tmpwA = np.ones_like(tmpA)
                    A.append(tmpA)
                    sA.append(tmpsA)
                    wA.append(abs(tmpwA))
                else:
                    print("Unexpected shape of tetramm data.")
    return A, sA, wA


def slitCalculator(energy=18.0,diodeSize=0.004,cameraLength=5.650):
    """
    Calculate slit width from given values...

    Energy in keV
    Measurements in m
    """
    twoTheta = np.arctan((diodeSize / 2) / cameraLength) * 2
    wavelength = (sp.h * sp.c) / (energy / 6.2415064799632e15)
    slitWidth = (4 * np.pi) / wavelength * np.sin(twoTheta / 2)
    # slitWidth = ((4 * np.pi) / ((sp.h * sp.c) / (wavelength))) * np.sin(twoTheta / 2)
    print("Using your specified diode size (" + str(diodeSize) + ") and camera length (" + str(cameraLength) + "),\nyour slit width is " + str(slitWidth * 1e-10) + " 1/A.")
    return slitWidth