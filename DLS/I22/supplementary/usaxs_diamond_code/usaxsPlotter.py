import os
import numpy as np
import matplotlib.pyplot as plt

plt.ioff()
from usaxsToolbox import directoryChecker, fileChecker, findNearest


def scatteringDataPlotter(scatteringDataObject, configDataObj, **kwargs):
    """
    Plot low and high gain scans of a scatteringDataObject.

    Various parameters can be specified including:
        xUnit: either "Q" (default) or "yaw"
        errBars: either None (default), "STD" or "SEM"
        xScale: "linear" (default) or "log"
        yScale: "linear" (default) or "log"
        xLims: None (default) or (lowerLimit, upperLimit)
        yLims: None (default) or (lowerLimit, upperLimit)
        legend: True (default) or a list of strings
        plotTitle: None (default) or a string
        fileName: None (default) or a string to an existing directory, with
                  file name appended
    If arguments which default to 'None' are not specified, these will not be
    used, i.e. no x limits will be set, no plot title will be set, and the
    figure will not be saved.
    """
    # Get various additional parameters (or lack thereof)
    xUnit = kwargs.get("xUnit", "Q")
    errBars = kwargs.get("errBars", False)
    mkrSize = kwargs.get("mkrSize", False)
    xScale = kwargs.get("xScale", "linear")
    yScale = kwargs.get("yScale", "linear")
    xLims = kwargs.get("xLims", None)
    yLims = kwargs.get("yLims", None)
    legend = kwargs.get("legend", True)
    plotTitle = kwargs.get("plotTitle", None)
    fileName = kwargs.get("fileName", None)

    # Setup
    xLo = {}
    xHi = {}

    yLo = {}
    yHi = {}

    if xUnit == "Q":
        for runNo in list(scatteringDataObject.usaxsRuns.keys()):
            xLo[runNo] = scatteringDataObject.usaxsRuns[runNo].loGain.Q
            xHi[runNo] = scatteringDataObject.usaxsRuns[runNo].hiGain.Q

            qMinX = [configDataObj.qMinTheo, configDataObj.qMinTheo]
            if scatteringDataObject.usaxsRuns[runNo].loGain.present:
                qVal1 = min([min(IData) for IData in scatteringDataObject.usaxsRuns[runNo].loGain.I])
                qVal2 = max([max(IData) for IData in scatteringDataObject.usaxsRuns[runNo].loGain.I])
            else:
                qVal1 = 0
                qVal2 = 1
            if scatteringDataObject.usaxsRuns[runNo].hiGain.present:
                qVal3 = min([min(IData) for IData in scatteringDataObject.usaxsRuns[runNo].hiGain.I])
                qVal4 = max([max(IData) for IData in scatteringDataObject.usaxsRuns[runNo].hiGain.I])
            else:
                qVal3 = 0
                qVal4 = 1
            qMinY = [np.min([qVal1, qVal3]), np.max([qVal2, qVal4])]

        xlabel = xUnit + "(" + configDataObj.outputUnit + ")"
    elif xUnit == "yaw":
        for runNo in list(scatteringDataObject.usaxsRuns.keys()):
            xLo[runNo] = scatteringDataObject.usaxsRuns[runNo].loGain.yaw
            xHi[runNo] = scatteringDataObject.usaxsRuns[runNo].hiGain.yaw

        xlabel = xUnit + "(" + configDataObj.inputUnit + ")"
    else:
        print("Unable to plot data, xUnit is not recognised.")

    if errBars is not False:
        # Extract errors. While errors for q do not exist before binning, errors
        # should always be present for I during fly scans (i.e. not a basic
        # TFG scan)

        xLoErr = {}
        xHiErr = {}
        yLoErr = {}
        yHiErr = {}

        if errBars == "STD":
            for runNo in list(scatteringDataObject.usaxsRuns.keys()):
                if any(
                    (
                        len(scatteringDataObject.usaxsRuns[runNo].loGain.QSTD) >= 1,
                        len(scatteringDataObject.usaxsRuns[runNo].hiGain.QSTD) >= 1,
                    )
                ):
                    xLoErr[runNo] = scatteringDataObject.usaxsRuns[runNo].loGain.QSTD
                    xHiErr[runNo] = scatteringDataObject.usaxsRuns[runNo].hiGain.QSTD
                else:
                    xLoErr[runNo] = np.full(len(xLo[runNo]), np.nan)
                    xHiErr[runNo] = np.full(len(xHi[runNo]), np.nan)
                if any(
                    (
                        len(scatteringDataObject.usaxsRuns[runNo].loGain.ISTD) >= 1,
                        len(scatteringDataObject.usaxsRuns[runNo].hiGain.ISTD) >= 1,
                    )
                ):
                    yLoErr[runNo] = [ISTD for ISTD in scatteringDataObject.usaxsRuns[runNo].loGain.ISTD]
                    yHiErr[runNo] = [ISTD for ISTD in scatteringDataObject.usaxsRuns[runNo].hiGain.ISTD]
                else:
                    print(
                        "Error bars cannot be plotted as no errors exist... Plotting normally."
                    )
                    errBars = False
        elif errBars == "SEM":
            for runNo in list(scatteringDataObject.usaxsRuns.keys()):
                if all(
                    (
                        len(scatteringDataObject.usaxsRuns[runNo].loGain.QSEM) >= 1,
                        len(scatteringDataObject.usaxsRuns[runNo].hiGain.QSEM) >= 1,
                    )
                ):
                    xLoErr[runNo] = scatteringDataObject.usaxsRuns[runNo].loGain.QSEM
                    xHiErr[runNo] = scatteringDataObject.usaxsRuns[runNo].hiGain.QSEM
                else:
                    xLoErr[runNo] = np.full(len(xLo[runNo]), np.nan)
                    xHiErr[runNo] = np.full(len(xHi[runNo]), np.nan)
                if all(
                    (
                        len(scatteringDataObject.usaxsRuns[runNo].loGain.ISEM) >= 1,
                        len(scatteringDataObject.usaxsRuns[runNo].hiGain.ISEM) >= 1,
                    )
                ):
                    yLoErr[runNo] = [ISEM for ISEM in scatteringDataObject.usaxsRuns[runNo].loGain.ISEM]
                    yHiErr[runNo] = [ISEM for ISEM in scatteringDataObject.usaxsRuns[runNo].hiGain.ISEM]
                else:
                    print(
                        "Error bars cannot be plotted as no errors exist... Plotting normally."
                    )
                    errBars = False
        else:
            print("Error bars name is not recognised... Plotting normally.")
            errBars = False

    yLo = {}
    yHi = {}

    for runNo in list(scatteringDataObject.usaxsRuns.keys()):
        yLo[runNo] = [IData for IData in scatteringDataObject.usaxsRuns[runNo].loGain.I]
        yHi[runNo] = [IData for IData in scatteringDataObject.usaxsRuns[runNo].hiGain.I]

    szLo = {}
    szHi = {}

    for runNo in list(scatteringDataObject.usaxsRuns.keys()):
        if mkrSize is not False:
            szLo[runNo] = [np.full_like(IData,mkrSize) for IData in yLo[runNo]]
            szHi[runNo] = [np.full_like(IData,mkrSize) for IData in yHi[runNo]]
        elif mkrSize is False:
            try:
                szLo[runNo] = [5 * (Iwt / max(Iwt)) for Iwt in scatteringDataObject.usaxsRuns[runNo].loGain.Iwt]  # normalise to 5  - a good size for viewing
            except:
                szLo[runNo] = [np.full_like(IData,5) for IData in yLo[runNo]]
            try:
                szHi[runNo] = [5 * (Iwt / max(Iwt)) for Iwt in scatteringDataObject.usaxsRuns[runNo].hiGain.Iwt] # normalaise to 5 - a good size for viewing
            except:
                szHi[runNo] = [np.full_like(IData,5) for IData in yHi[runNo]]

    # Plot
    plt.figure(figsize=(8, 6), dpi=172)

    if not errBars:
        for runNo in list(scatteringDataObject.usaxsRuns.keys()):
            for inc in range(len(yLo[runNo])):
                plt.scatter(
                    xLo[runNo],
                    yLo[runNo][inc],
                    s=szLo[runNo][inc],
                    alpha=0.5,
                    label=(f"low gain: Run{runNo}, Diode{inc}"),
                )
            for inc in range(len(yHi[runNo])):
                plt.scatter(
                    xHi[runNo],
                    yHi[runNo][inc],
                    s=szHi[runNo][inc],
                    alpha=0.5,
                    label=(f"high gain: Run{runNo}, Diode{inc}"),
                )
    else:
        for runNo in list(scatteringDataObject.usaxsRuns.keys()):
            for inc in range(len(yLo[runNo])):
                plt.errorbar(
                    xLo[runNo],
                    yLo[runNo][inc],
                    xerr=xLoErr[runNo],
                    yerr=yLoErr[runNo][inc],
                    ms=3,
                    fmt="o",
                    alpha=0.5,
                    label=(f"low gain: Run{runNo}, Diode{inc}"),
                )
            for inc in range(len(yHi[runNo])):
                plt.errorbar(
                    xHi[runNo],
                    yHi[runNo][inc],
                    xerr=xHiErr[runNo],
                    yerr=yHiErr[runNo][inc],
                    ms=3,
                    fmt="o",
                    alpha=0.5,
                    label=(f"low gain: Run{runNo}, Diode{inc}"),
            )

    if xUnit == "Q":
        plt.plot(
            qMinX,
            qMinY,
            color="#232323",
            alpha=0.5,
            label="theoretical q min",
        )

    plt.ylabel("Intensity")
    plt.xlabel(xlabel)
    plt.xscale(xScale)
    plt.yscale(yScale)

    if xLims is not None:
        plt.xlim(xLims)

    if yLims is not None:
        plt.ylim(yLims)

    if plotTitle is not None:
        plt.title(plotTitle)

    if legend is not False:
        if legend is True:
            plt.legend()
        else:
            plt.legend(legend)

    if fileName is not None:
        dirName = "/".join(fileName.split("/")[:-1])
        directoryChecker(dirName)
        fileChecker(fileName)
        plt.savefig(fileName)

    # plt.show()
    # plt.close()


def foldedDataPlotter(scatteringDataObject, configDataObj, **kwargs):
    """
    Plot folded low and high gain scans of a scatteringDataObject.

    Various parameters can be specified including:
        xUnit: either "Q" (default) or "yaw"
        errBars: either None (default), "STD" or "SEM"
        xScale: "linear" (default) or "log"
        yScale: "linear" (default) or "log"
        xLims: None (default) or (lowerLimit, upperLimit)
        yLims: None (default) or (lowerLimit, upperLimit)
        legend: True (default) or a list of strings
        plotTitle: None (default) or a string
        fileName: None (default) or a string to an existing directory, with
                  file name appended
    If arguments which default to 'None' are not specified, these will not be
    used, i.e. no x limits will be set, no plot title will be set, and the
    figure will not be saved.
    """
    # Get various additional parameters (or lack thereof)
    xUnit = kwargs.get("xUnit", "Q")
    mkrSize = kwargs.get("mkrSize", False)
    xScale = kwargs.get("xScale", "linear")
    yScale = kwargs.get("yScale", "linear")
    legend = kwargs.get("legend", True)
    plotTitle = kwargs.get("plotTitle", None)
    fileName = kwargs.get("fileName", None)

    # Setup
    xLo = {}
    xHi = {}

    yLo = {}
    yHi = {}

    if xUnit == "Q":
        for runNo in list(scatteringDataObject.usaxsRuns.keys()):
            xLo[runNo] = scatteringDataObject.usaxsRuns[runNo].loGain.Q
            xHi[runNo] = scatteringDataObject.usaxsRuns[runNo].hiGain.Q
            
            qMinX = [configDataObj.qMinTheo, configDataObj.qMinTheo]
            if scatteringDataObject.usaxsRuns[runNo].loGain.present:
                qVal1 = min([min(IData) for IData in scatteringDataObject.usaxsRuns[runNo].loGain.I])
                qVal2 = max([max(IData) for IData in scatteringDataObject.usaxsRuns[runNo].loGain.I])
            else:
                qVal1 = 0
                qVal2 = 1
            if scatteringDataObject.usaxsRuns[runNo].hiGain.present:
                qVal3 = min([min(IData) for IData in scatteringDataObject.usaxsRuns[runNo].hiGain.I])
                qVal4 = max([max(IData) for IData in scatteringDataObject.usaxsRuns[runNo].hiGain.I])
            else:
                qVal3 = 0
                qVal4 = 1
            qMinY = [np.min([qVal1, qVal3]), np.max([qVal2, qVal4])]

        xlabel = xUnit + "(" + configDataObj.outputUnit + ")"
    elif xUnit == "yaw":
        for runNo in list(scatteringDataObject.usaxsRuns.keys()):
            xLo[runNo] = scatteringDataObject.usaxsRuns[runNo].loGain.yaw
            xHi[runNo] = scatteringDataObject.usaxsRuns[runNo].hiGain.yaw

        xlabel = xUnit + "(" + configDataObj.inputUnit + ")"
    else:
        print("Unable to plot data, xUnit is not recognised.")

    yLo = {}
    yHi = {}

    for runNo in list(scatteringDataObject.usaxsRuns.keys()):
        yLo[runNo] = [IData for IData in scatteringDataObject.usaxsRuns[runNo].loGain.I]
        yHi[runNo] = [IData for IData in scatteringDataObject.usaxsRuns[runNo].hiGain.I]

    szLo = {}
    szHi = {}

    for runNo in list(scatteringDataObject.usaxsRuns.keys()):
        if mkrSize is not False:
            szLo[runNo] = mkrSize
            szHi[runNo] = mkrSize
        elif mkrSize is False:
            try:
                szLo[runNo] = 5 * (
                    scatteringDataObject.usaxsRuns[runNo].loGain.Iwt[0]
                    / max(scatteringDataObject.usaxsRuns[runNo].loGain.Iwt[0])
                )  # normalise to 5  - a good size for viewing
            except:
                szLo[runNo] = 5
            try:
                szHi[runNo] = 5 * (
                    scatteringDataObject.usaxsRuns[runNo].hiGain.Iwt[0]
                    / max(scatteringDataObject.usaxsRuns[runNo].hiGain.Iwt[0])
                )  # normalaise to 5 - a good size for viewing
            except:
                szHi[runNo] = 5
        # elif mkrSize is False:
        #     szLo[runNo] = 5 * (scatteringDataObject.usaxsRuns[runNo].loGain.Iwt / max(scatteringDataObject.usaxsRuns[runNo].loGain.Iwt)) # normalise to 5  - a good size for viewing
        #     szHi[runNo] = 5 * (scatteringDataObject.usaxsRuns[runNo].hiGain.Iwt / max(scatteringDataObject.usaxsRuns[runNo].hiGain.Iwt)) # normalaise to 5 - a good size for viewing

    # Plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 6), dpi=172, sharey=True)

    for runNo in list(scatteringDataObject.usaxsRuns.keys()):
        if scatteringDataObject.usaxsRuns[runNo].loGain.present:
            for inc in range(len(yLo[runNo])):
                ax1.scatter(
                    abs(xLo[runNo][xLo[runNo] < 0]),
                    yLo[runNo][inc][xLo[runNo] < 0],
                    s=szLo[runNo][xLo[runNo] < 0],
                    c="#440154",
                    alpha=0.5,
                    label=(f"- low gain: Run{runNo}, Diode{inc}"),
                )
                ax1.scatter(
                    abs(xLo[runNo][xLo[runNo] > 0]),
                    yLo[runNo][inc][xLo[runNo] > 0],
                    s=szLo[runNo][xLo[runNo] > 0],
                    c="#73D055",
                    alpha=0.5,
                    label=(f"+ low gain: Run{runNo}, Diode{inc}"),
                )
                ax2.scatter(
                    abs(xLo[runNo][xLo[runNo] < 0]),
                    yLo[runNo][inc][xLo[runNo] < 0],
                    s=szLo[runNo][xLo[runNo] < 0],
                    c="#440154",
                    alpha=0.5,
                    label=(f"- low gain: Run{runNo}, Diode{inc}"),
                )
                ax2.scatter(
                    abs(xLo[runNo][xLo[runNo] > 0]),
                    yLo[runNo][inc][xLo[runNo] > 0],
                    s=szLo[runNo][xLo[runNo] > 0],
                    c="#73D055",
                    alpha=0.5,
                    label=(f"+ low gain: Run{runNo}, Diode{inc}"),
                )
        for inc in range(len(yHi[runNo])):
            if scatteringDataObject.usaxsRuns[runNo].hiGain.present:
                ax1.scatter(
                    abs(xHi[runNo][xHi[runNo] < 0]),
                    yHi[runNo][inc][xHi[runNo] < 0],
                    s=szHi[runNo][xHi[runNo] < 0],
                    c="#33638D",
                    alpha=0.5,
                    label=(f"- high gain: Run{runNo}, Diode{inc}"),
                )
                ax1.scatter(
                    abs(xHi[runNo][xHi[runNo] > 0]),
                    yHi[runNo][inc][xHi[runNo] > 0],
                    s=szHi[runNo][xHi[runNo] > 0],
                    c="#DCE319",
                    alpha=0.5,
                    label=(f"+ high gain: Run{runNo}, Diode{inc}"),
                )
                ax2.scatter(
                    abs(xHi[runNo][xHi[runNo] < 0]),
                    yHi[runNo][inc][xHi[runNo] < 0],
                    s=szHi[runNo][xHi[runNo] < 0],
                    c="#33638D",
                    alpha=0.5,
                    label=(f"- high gain: Run{runNo}, Diode{inc}"),
                )
                ax2.scatter(
                    abs(xHi[runNo][xHi[runNo] > 0]),
                    yHi[runNo][inc][xHi[runNo] > 0],
                    s=szHi[runNo][xHi[runNo] > 0],
                    c="#DCE319",
                    alpha=0.5,
                    label=(f"+ high gain: Run{runNo}, Diode{inc}"),
                )

    if xUnit == "Q":
        plt.plot(
            qMinX,
            qMinY,
            color="#232323",
            alpha=0.5,
            label="theoretical q min",
        )

    ax1.set_ylabel("Intensity")
    ax1.set_xlabel(xlabel)
    ax2.set_xlabel(xlabel)
    ax1.set_xscale("linear")
    ax2.set_xscale("log")
    ax2.set_xlim((1e0, 1e2))
    ax1.set_yscale(yScale)
    ax1.set_title("Full Range")
    ax2.set_title("Oxidation Check")

    if plotTitle is not None:
        fig.suptitle(plotTitle)

    if legend is not False:
        if legend is True:
            ax2.legend()
        else:
            ax2.legend(legend)

    if fileName is not None:
        dirName = "/".join(fileName.split("/")[:-1])
        directoryChecker(dirName)
        fileChecker(fileName)
        plt.savefig(fileName)

    # plt.show()
    # plt.close()


def mergedDataPlotter(scatteringDataObject, configDataObj, **kwargs):
    """
    Plot merged data of a scatteringDataObject.

    Various parameters can be specified including:
        xUnit: either "Q" (default) or "yaw"
        errBars: either None (default), "STD" or "SEM"
        xScale: "linear" (default) or "log"
        yScale: "linear" (default) or "log"
        xLims: None (default) or (lowerLimit, upperLimit)
        yLims: None (default) or (lowerLimit, upperLimit)
        legend: True (default) or a list of strings
        plotTitle: None (default) or a string
        fileName: None (default) or a string to an existing directory, with
                  file name appended
    If arguments which default to 'None' are not specified, these will not be
    used, i.e. no x limits will be set, no plot title will be set, and the
    figure will not be saved.
    """
    # Get various additional parameters (or lack thereof)
    xUnit = kwargs.get("xUnit", "Q")
    errBars = kwargs.get("errBars", False)
    xScale = kwargs.get("xScale", "linear")
    yScale = kwargs.get("yScale", "linear")
    xLims = kwargs.get("xLims", None)
    yLims = kwargs.get("yLims", None)
    legend = kwargs.get("legend", True)
    plotTitle = kwargs.get("plotTitle", None)
    fileName = kwargs.get("fileName", None)

    # Setup
    if xUnit == "Q":
        x = scatteringDataObject.mData.Q

        qMinX = [configDataObj.qMinTheo, configDataObj.qMinTheo]
        qMinY = [min(scatteringDataObject.mData.I), max(scatteringDataObject.mData.I)]

        xlabel = xUnit + "(" + configDataObj.outputUnit + ")"
    elif xUnit == "yaw":
        x = scatteringDataObject.mData.yaw

        xlabel = xUnit + "(" + configDataObj.inputUnit + ")"
    else:
        print("Unable to plot data, xUnit is not recognised.")

    if errBars is not False:
        if (
            len(scatteringDataObject.mData.QSTD) >= 1
            and len(scatteringDataObject.mData.ISTD) >= 1
        ):
            if errBars == "STD":
                xErr = scatteringDataObject.mData.QSTD
                yErr = scatteringDataObject.mData.ISTD
            elif errBars == "SEM":
                xErr = scatteringDataObject.mData.QSEM
                yErr = scatteringDataObject.mData.ISEM
            else:
                print("Error bars name is not recognised... Plotting normally.")
                errBars = False
        else:
            print(
                "Error bars cannot be plotted as no errors exist... Plotting normally."
            )
            errBars = False

    y = scatteringDataObject.mData.I

    # Plot
    plt.figure(figsize=(8, 6), dpi=172)

    if not errBars:
        plt.scatter(x, y, s=3, alpha=0.5, label="merged data")
    else:
        plt.errorbar(
            x, y, xerr=xErr, yerr=yErr, ms=3, fmt="o", alpha=0.5, label="merged data"
        )

    if xUnit == "Q":
        plt.plot(
            qMinX,
            qMinY,
            color="#232323",
            alpha=0.5,
            label="theoretical q min",
        )

    plt.ylabel("Intensity")
    plt.xlabel(xlabel)
    plt.xscale(xScale)
    plt.yscale(yScale)

    if xLims is not None:
        plt.xlim(xLims)

    if yLims is not None:
        plt.ylim(yLims)

    if plotTitle is not None:
        plt.title(plotTitle)

    if legend is not False:
        if legend is True:
            plt.legend()
        else:
            plt.legend(legend)

    if fileName is not None:
        dirName = "/".join(fileName.split("/")[:-1])
        directoryChecker(dirName)
        fileChecker(fileName)
        plt.savefig(fileName)

    # plt.show()
    # plt.close()


def saxsDataPlotter(scatteringDataObject, configDataObj, **kwargs):
    """
    Plot merged data of a scatteringDataObject.

    Various parameters can be specified including:
        xUnit: either "Q" (default) or "yaw"
        errBars: either None (default), "STD" or "SEM"
        xScale: "linear" (default) or "log"
        yScale: "linear" (default) or "log"
        xLims: None (default) or (lowerLimit, upperLimit)
        yLims: None (default) or (lowerLimit, upperLimit)
        legend: True (default) or a list of strings
        plotTitle: None (default) or a string
        fileName: None (default) or a string to an existing directory, with
                  file name appended
    If arguments which default to 'None' are not specified, these will not be
    used, i.e. no x limits will be set, no plot title will be set, and the
    figure will not be saved.
    """
    # Get various additional parameters (or lack thereof)
    xUnit = kwargs.get("xUnit", "Q")
    errBars = kwargs.get("errBars", False)
    xScale = kwargs.get("xScale", "linear")
    yScale = kwargs.get("yScale", "linear")
    xLims = kwargs.get("xLims", None)
    yLims = kwargs.get("yLims", None)
    legend = kwargs.get("legend", True)
    plotTitle = kwargs.get("plotTitle", None)
    fileName = kwargs.get("fileName", None)

    # Setup
    if xUnit == "Q":
        x = scatteringDataObject.mData.Q

        qMinX = [configDataObj.qMinTheo]
        idx1 = findNearest(scatteringDataObject.mData.Q, qMinX)
        qMinY = scatteringDataObject.mData.I[idx1]

        xlabel = xUnit + "(" + configDataObj.outputUnit + ")"
    elif xUnit == "yaw":
        print("Unable to plot merged data against yaw values.")
    else:
        print("Unable to plot data, xUnit is not recognised.")

    if errBars is not False:
        if (
            len(scatteringDataObject.mData.QSTD) >= 1
            and len(scatteringDataObject.mData.ISTD) >= 1
        ):
            if errBars == "STD":
                xErr = scatteringDataObject.mData.QSTD
                yErr = scatteringDataObject.mData.ISTD
            elif errBars == "SEM":
                xErr = scatteringDataObject.mData.QSEM
                yErr = scatteringDataObject.mData.ISEM
            else:
                print("Error bars name is not recognised... Plotting normally.")
                errBars = False
        else:
            print(
                "Error bars cannot be plotted as no errors exist... Plotting normally."
            )
            errBars = False

    y = scatteringDataObject.mData.I

    if scatteringDataObject.mData.desmearedQ != []:
        plt_dsm = True
        x_dsm = scatteringDataObject.mData.desmearedQ
        y_dsm = scatteringDataObject.mData.desmearedI
        yErr_dsm = scatteringDataObject.mData.desmearedIerr
    else:
        plt_dsm = False

    saxsx = scatteringDataObject.mData.saxsQ
    saxsy = scatteringDataObject.mData.saxsI
    saxsyerr = scatteringDataObject.mData.saxsIerr

    # Plot
    plt.figure(figsize=(8, 6), dpi=172)

    if errBars is False:
        plt.scatter(x, y, s=3, alpha=0.5, label="USAXS data")
        plt.scatter(saxsx, saxsy, s=3, alpha=0.5, label="SAXS data")
        if plt_dsm:
            plt.scatter(x_dsm, y_dsm, s=3, alpha=0.5, label="desmeared USAXS data")
    else:
        plt.errorbar(
            x, y, xerr=xErr, yerr=yErr, ms=3, fmt="o", alpha=0.5, label="USAXS data"
        )
        plt.errorbar(
            saxsx, saxsy, yerr=saxsyerr, ms=3, fmt="o", alpha=0.5, label="SAXS data"
        )
        if plt_dsm:
            plt.errorbar(
                x_dsm,
                y_dsm,
                yerr=yErr_dsm,
                ms=3,
                fmt="o",
                alpha=0.5,
                label="desmeared USAXS data",
            )

    if xUnit == "Q":
        plt.scatter(
            qMinX,
            qMinY,
            s=30,
            color="#232323",
            marker="o",
            alpha=1,
            label="theoretical q min",
        )

    plt.ylabel("Intensity")
    plt.xlabel(xlabel)
    plt.xscale(xScale)
    plt.yscale(yScale)

    if xLims is not None:
        plt.xlim(xLims)

    if yLims is not None:
        plt.ylim(yLims)

    if plotTitle is not None:
        plt.title(plotTitle)

    if legend is not False:
        if legend is True:
            plt.legend()
        else:
            plt.legend(legend)

    if fileName is not None:
        dirName = "/".join(fileName.split("/")[:-1])
        directoryChecker(dirName)
        fileChecker(fileName)
        plt.savefig(fileName)

    # plt.show()
    # plt.close()
