import glob
from nexusformat.nexus import nxopen, NXdata, NXentry, NXfield, NXlink
import numpy as np
from usaxsDataClasses import scatteringDataObj, rawUsaxsDataObj
from usaxsToolbox import directoryChecker, fileChecker


def scanNosGrabber(scanNumber, f, feedback):
    try:
        if feedback:
            print("Trying to find scan numbers at entry.sample.scanNos")
        scanNosStr = str(f["entry/sample/scanNos"])[1:-1]
        scanNos = scanNosStr.split(", ")
        scanNos.append(str(scanNumber))  # add this scan as well - it's not in the list!
    except:
        scanNos = [str(scanNumber)]
    if feedback:
        print(f"Associated scan numbers in list: {scanNos}")
    return scanNos


def divvyScanNos(scanNos, configDataObj):
    """Divvy up scan numbers based on metadata."""
    feedback = configDataObj.feedback
    dataPath = configDataObj.dataPath
    scanNosAll = []
    scanNosUnique = []

    for scanNo in scanNos:
        if feedback:
            print(f"\nRetrieving metadata for scan number {scanNo}")
        fileName = (
            dataPath
            + "/"
            + configDataObj.fileNamePrefix
            + str(scanNo)
            + configDataObj.fileNameSuffix
        )
        try:
            with nxopen(fileName) as f:
                scanNosAll.append(scanNosGrabber(scanNo, f, feedback))
                if feedback:
                    print(f"Retrieved metadata for scan number {scanNo}")
        except:
            print(f"!! Could not retrieve metadata for scan number {scanNo}")

    scanNosAll.sort(
        key=len, reverse=True
    )  # sort by length of list so longest are first
    for someScanNos in scanNosAll:
        is_subset = False
        for otherScanNos in scanNosAll:
            if someScanNos != otherScanNos and set(someScanNos).issubset(
                set(otherScanNos)
            ):
                is_subset = True
                break
        if not is_subset:
            scanNosUnique.append(someScanNos)
    return scanNosUnique


def nexusNodeGrabber(f, nexusPaths, feedback):
    """Retrieve some data from a set of nexus paths.

    Args:
        f (NXroot): NXroot object containing the NeXus tree
        nexusPaths (list): list of lists containing (potential) paths to requested datasets
        feedback (bool): bool turning verbose feedback on/off
    """
    data = []
    for datasetNexusPaths in nexusPaths:
        dat = None
        for nexusPath in datasetNexusPaths:
            if (
                not dat
            ):  # As soon as you have the data, you don't have to look for it anymore
                try:
                    if feedback:
                        print(f"Trying to find data at {nexusPath}...")
                    dat = f[nexusPath]
                    if feedback:
                        print(f"Found data at {nexusPath}.")
                except:
                    if feedback:
                        print(f"!! Failed to find data at {nexusPath}...")
        data.append(dat)
    return data


def scanTitleGrabber(f, nexusPaths, feedback):
    data = nexusNodeGrabber(f, nexusPaths, feedback)
    sampleName = [str(name) for name in data][0]  # You will only ever have one
    return sampleName


def yawDataGrabber(f, nexusPaths, feedback):
    data = nexusNodeGrabber(f, nexusPaths, feedback)
    yawData = [np.squeeze(np.array(dat)) for dat in data][
        0
    ]  # You will only ever have one
    return yawData


def tetrammDataGrabber(f, nexusPaths, channels, feedback):
    data = nexusNodeGrabber(f, nexusPaths, feedback)
    diodeData = []
    for dat, channel in zip(data, channels):
        diodeData.append(np.squeeze(np.array(dat))[:, :, channel])
    return diodeData


def countTimeGrabber(f, nexusPaths, feedback):
    try:
        data = nexusNodeGrabber(f, nexusPaths, feedback)
    except:
        data = [1.0]
        if feedback:
            print("!! Count time not found, defaulting to 1.0")
    countTime = [np.squeeze(np.array(dat)) for dat in data][
        0
    ]  # You will only ever have one
    return countTime


def processedSAXSDataGrabber(processedSAXSDataPath, scanNumber):
    """Retrieve q, I, Ierr from pre-processed SAXS data.

    Attempt to retrieve data from pre-processed SAXS data, from provided path
    and known scan number. If > 1 data file exists for the scan number, the
    most recent is chosen. Any frames are averaged.
    """

    def saxsDataGrabber(dataPath):
        with nxopen(dataPath) as fSAXS:
            q = np.array(fSAXS["processed/result/q"])
            data = np.squeeze(np.squeeze(np.array(fSAXS["processed/result/data"])))
            if data.ndim == 2:
                data = np.average(data, axis=0)
            errors = np.squeeze(np.squeeze(np.array(fSAXS["processed/result/errors"])))
            if errors.ndim == 2:
                errors = np.average(errors, axis=0)
        return q, data, errors

    fileNames = glob.glob((processedSAXSDataPath + "/*" + str(scanNumber) + "*.nxs"))
    if len(fileNames) == 1:
        print(
            f"Processed SAXS data for scan number {scanNumber} was found, processing now..."
        )
        dataPath = fileNames[0]
    elif len(fileNames) > 1:
        # check if glob.glob((processedSAXSDataPath + "/*" + str(scanNumber) + "*saxs*.nxs")) is only one
        # If there's still more than one, then get the most recent one
        print(
            f"More than one file for processed SAXS data for scan number {scanNumber} was found, processing now the most recent file..."
        )
        timeDate = np.array([])
        inc = 0
        for file in fileNames:
            timeDate[inc] = int("".join(file.split("_")[-2:-1]))
            inc += 1
        maxArg = np.maxarg(timeDate)
        dataPath = fileNames[maxArg]
    q, data, errors = saxsDataGrabber(dataPath)
    return q, data, errors


def genericDataGrabber(scanNos, dataPath, configDataObj, excludeString, includeString):
    """Retrieves selected data from a list of scan numbers.

    Attempts to retrieve datasets from I0, PANDA_05, BSDIODES, and
    sample name, for a selected scan number. If not all datasets (aside from I0)
    can be retrieved, no datasets for the selected scan number
    are saved, and user is notified. Returns a dictionary of scatteringDataObj,
    with each unique sample name making up each scatteringDataObj.
    An 'excludeString' and an 'includeString' can be set, to triage scans by
    sample name. If 'excludeString' is found, scans are discarded. If
    'includeString' is not found, scans are discarded. *Useful for finding dark
    current datasets in a large number of scans.
    """
    feedback = configDataObj.feedback
    scatteringData = {}

    for scanNo in scanNos:
        print(f"\nRetrieving data for scan number {scanNo}.")
        fileName = (
            dataPath
            + "/"
            + configDataObj.fileNamePrefix
            + str(scanNo)
            + configDataObj.fileNameSuffix
        )
        try:
            with nxopen(fileName) as f:
                try:
                    nexusPaths = configDataObj.titleNexusPath
                    scanTitle = scanTitleGrabber(f, nexusPaths, feedback)

                    # Check for inclusions / exclusions with your strings
                    if excludeString is None:
                        if includeString is None:
                            getScan = True  # No exclusion and no inclusion
                        elif includeString is not None:
                            if includeString in scanTitle:
                                getScan = True  # No exclusion and inclusion found
                            elif includeString not in scanTitle:
                                getScan = False  # No exclusion and inclusion not found
                    elif excludeString is not None:
                        if excludeString in scanTitle:
                            getScan = False  # Exclusion found
                        if excludeString not in scanTitle:
                            if includeString is None:
                                getScan = True  # Exclusion not found and no inclusion
                            elif includeString is not None:
                                if includeString in scanTitle:
                                    getScan = (
                                        True  # Exclusion not found and inclusion found
                                    )
                                elif includeString not in scanTitle:
                                    getScan = False  # Exclusion not found and inclusion not found

                    if getScan:
                        if " usaxs scan" in scanTitle:
                            scanType = "usaxs"
                            sampleName = scanTitle.removesuffix(" usaxs scan")
                            if feedback:
                                print(f"Scan number {scanNo} is a USAXS scan.")

                        elif " dc scan" in scanTitle:
                            scanType = "dc"
                            sampleName = scanTitle.removesuffix(" dc scan")
                            if feedback:
                                print(f"Scan number {scanNo} is a dark current scan.")

                        elif "SWAXS scan" in scanTitle:
                            scanType = "swaxs"
                            sampleName = scanTitle.removesuffix(" SWAXS scan")
                            runNo = 0
                            if feedback:
                                print(f"Scan number {scanNo} is a SWAXS scan.")

                        else:
                            scanType = "unknown"
                            sampleName = scanTitle
                            runNo = 0
                            if feedback:
                                print(f"Scan number {scanNo} is of unknown scan type.")

                        if "low gain" in sampleName:
                            scanGain = "low"
                            sampleName = sampleName.removesuffix("low gain")
                            if feedback:
                                print(f"Scan number {scanNo} is a low gain scan.")

                        elif "high gain" in sampleName:
                            scanGain = "high"
                            sampleName = sampleName.removesuffix("high gain")
                            if feedback:
                                print(f"Scan number {scanNo} is a high gain scan.")

                        sampleName = sampleName.strip(
                            " -:"
                        )  # Remove unwanted trailing characters which may be there...

                        try:
                            if any((scanType == "usaxs", scanType == "dc")):
                                if feedback:
                                    print("Retrieving run number from scan title")
                                runNo = sampleName.split("Run: ")[1]
                                sampleName = sampleName.split("Run: ")[0]
                                if feedback:
                                    print(f"Usaxs scan is run number {runNo}")
                        except:
                            if feedback:
                                print(
                                    "!! Run number could not be determined, defaulting to '0'"
                                )
                            runNo = "0"

                        sampleName = sampleName.strip(
                            " -:"
                        )  # Remove unwanted trailing characters which may be there...

                        if scanType != "unknown":
                            try:
                                scatteringData[sampleName]
                                if feedback:
                                    print(
                                        f"scatteringDataObject for {sampleName} already exists..."
                                    )
                            except:
                                scatteringData[sampleName] = scatteringDataObj()
                                scatteringData[sampleName].setSampleName(sampleName)
                                if feedback:
                                    print(
                                        f"scatteringDataObject for {sampleName} does not exist, creating now..."
                                    )

                            if any((scanType == "usaxs", scanType == "dc")):
                                if feedback:
                                    print(
                                        f"Scan number {scanNo} is a dark current or USAXS scan, checking for rawUsaxsDataObj now..."
                                    )
                                try:
                                    scatteringData[sampleName].usaxsRuns[runNo]
                                    if feedback:
                                        print(
                                            f"rawUsaxsDataObj for run number {runNo} of {sampleName} already exists..."
                                        )
                                except:
                                    scatteringData[sampleName].usaxsRuns[runNo] = (
                                        rawUsaxsDataObj()
                                    )
                                    scatteringData[sampleName].usaxsRuns[
                                        runNo
                                    ].setRunNo(runNo)
                                    if feedback:
                                        print(
                                            f"rawUsaxsDataObj for run number {runNo} of {sampleName} does not exist, creating now..."
                                        )

                            # Read in dark current data
                            if scanType == "dc":
                                print(
                                    f"Scan number {scanNo} is a dark current scan, grabbing data now..."
                                )
                                nexusPaths = configDataObj.diodeNexusPath
                                channels = configDataObj.diodeChannels
                                diodeData = tetrammDataGrabber(
                                    f=f,
                                    nexusPaths=nexusPaths,
                                    channels=channels,
                                    feedback=feedback,
                                )
                                nexusPaths = configDataObj.countNexusPath
                                countTime = countTimeGrabber(
                                    f=f, nexusPaths=nexusPaths, feedback=feedback
                                )

                                if scanGain == "low":
                                    scatteringData[sampleName].usaxsRuns[
                                        runNo
                                    ].loGainDC.countTime = countTime
                                    scatteringData[sampleName].usaxsRuns[
                                        runNo
                                    ].loGainDC.setIRaw(diodeData)
                                    scatteringData[sampleName].usaxsRuns[
                                        runNo
                                    ].loGainDC.setScanNo(scanNo)
                                    scatteringData[sampleName].usaxsRuns[
                                        runNo
                                    ].loGainDC.present = True
                                    if feedback:
                                        print(
                                            "Low gain dark current scan data acquired"
                                        )
                                elif scanGain == "high":
                                    scatteringData[sampleName].usaxsRuns[
                                        runNo
                                    ].hiGainDC.countTime = countTime
                                    scatteringData[sampleName].usaxsRuns[
                                        runNo
                                    ].hiGainDC.setIRaw(diodeData)
                                    scatteringData[sampleName].usaxsRuns[
                                        runNo
                                    ].hiGainDC.setScanNo(scanNo)
                                    scatteringData[sampleName].usaxsRuns[
                                        runNo
                                    ].hiGainDC.present = True
                                    if feedback:
                                        print(
                                            "High gain dark current scan data acquired"
                                        )

                            # Read in USAXS data
                            elif scanType == "usaxs":
                                print(
                                    f"Scan number {scanNo} is a USAXS scan, grabbing data now..."
                                )
                                nexusPaths = configDataObj.yawNexusPath
                                yawData = yawDataGrabber(
                                    f=f, nexusPaths=nexusPaths, feedback=feedback
                                )
                                nexusPaths = configDataObj.diodeNexusPath
                                channels = configDataObj.diodeChannels
                                diodeData = tetrammDataGrabber(
                                    f=f,
                                    nexusPaths=nexusPaths,
                                    channels=channels,
                                    feedback=feedback,
                                )
                                nexusPaths = configDataObj.i0NexusPath
                                channels = configDataObj.i0Channels
                                i0Data = tetrammDataGrabber(
                                    f=f,
                                    nexusPaths=nexusPaths,
                                    channels=channels,
                                    feedback=feedback,
                                )
                                nexusPaths = configDataObj.countNexusPath
                                countTime = countTimeGrabber(
                                    f=f, nexusPaths=nexusPaths, feedback=feedback
                                )

                                if scanGain == "low":
                                    scatteringData[sampleName].usaxsRuns[
                                        runNo
                                    ].loGain.countTime = countTime
                                    scatteringData[sampleName].usaxsRuns[
                                        runNo
                                    ].loGain.setYawRaw(yawData)
                                    scatteringData[sampleName].usaxsRuns[
                                        runNo
                                    ].loGain.setIRaw(diodeData)
                                    scatteringData[sampleName].usaxsRuns[
                                        runNo
                                    ].loGain.setI0Raw(i0Data)
                                    scatteringData[sampleName].usaxsRuns[
                                        runNo
                                    ].loGain.setScanNo(scanNo)
                                    scatteringData[sampleName].usaxsRuns[
                                        runNo
                                    ].loGain.present = True
                                    if feedback:
                                        print("Low gain usaxs scan data acquired")
                                elif scanGain == "high":
                                    scatteringData[sampleName].usaxsRuns[
                                        runNo
                                    ].hiGain.countTime = countTime
                                    scatteringData[sampleName].usaxsRuns[
                                        runNo
                                    ].hiGain.setYawRaw(yawData)
                                    scatteringData[sampleName].usaxsRuns[
                                        runNo
                                    ].hiGain.setIRaw(diodeData)
                                    scatteringData[sampleName].usaxsRuns[
                                        runNo
                                    ].hiGain.setI0Raw(i0Data)
                                    scatteringData[sampleName].usaxsRuns[
                                        runNo
                                    ].hiGain.setScanNo(scanNo)
                                    scatteringData[sampleName].usaxsRuns[
                                        runNo
                                    ].hiGain.present = True
                                    if feedback:
                                        print("High gain usaxs scan data acquired")
                                scatteringData[sampleName].usaxsPresent = True

                            # Read in SWAXS data
                            elif scanType == "swaxs":
                                print(
                                    f"Scan number {scanNo} is a SAXS-WAXS frame, searching for processed SAXS data..."
                                )
                                try:
                                    q, data, errors = processedSAXSDataGrabber(
                                        configDataObj.processedSAXSDataPath, scanNo
                                    )
                                    scatteringData[sampleName].mData.saxsScanNo = scanNo
                                    scatteringData[sampleName].mData.saxsI = data
                                    scatteringData[sampleName].mData.saxsIerr = errors
                                    scatteringData[sampleName].mData.saxsQ = q
                                    scatteringData[sampleName].swaxsPresent = True
                                except:
                                    print(
                                        f"!! Could not retrieve processed SAXS data for scan number {scanNo} moving on..."
                                    )

                        elif scanType == "unknown":
                            print(
                                f"Scan type for scan number {scanNo} is unknown (i.e. not a dark current, USAXS, or SAXS/WAXS scan, moving on..."
                            )

                    else:
                        print(
                            f"Scan number {scanNo} excluded because either excludeString ('{excludeString}') was found or includeString ('{includeString}') was not found..."
                        )
                except:
                    print(
                        f"!! Failed to retrieve necessary data for scan number {scanNo}, moving on..."
                    )
        except:
            print(f"!! Data from scan number {scanNo} could not be retrieved.")

    return scatteringData


def backgroundDataGrabber(configDataObj, **kwargs):
    """Retrieves background data as defined in configDataObj.

    An 'excludeString' and an 'includeString' can be set, to triage scans by
    sample name. If 'excludeString' is found, scans are discarded. If
    'includeString' is not found, scans are discarded.
    """
    excludeString = kwargs.get("excludeString", None)
    includeString = kwargs.get("includeString", None)

    print("Now retrieving background data...")

    scanNos = configDataObj.backgroundNos
    if len(scanNos) >= 1:
        dataPath = configDataObj.dataPath
        return genericDataGrabber(
            scanNos,
            dataPath,
            configDataObj,
            excludeString=excludeString,
            includeString=includeString,
        )
    else:
        print("No file numbers found to retrieve background data")
        pass


def sampleDataGrabber(configDataObj, **kwargs):
    """Retrieves sample data as defined in configDataObj.

    An 'excludeString' and an 'includeString' can be set, to triage scans by
    sample name. If 'excludeString' is found, scans are discarded. If
    'includeString' is not found, scans are discarded. *Useful for disregarding
    dark current scans when reading in sample data (i.e. if dark current scans
    are throughout your samples).
    """
    excludeString = kwargs.get("excludeString", None)
    includeString = kwargs.get("includeString", None)

    print("Now retrieving sample data...")

    scanNos = configDataObj.scanNos
    if len(scanNos) >= 1:
        dataPath = configDataObj.dataPath
        return genericDataGrabber(
            scanNos,
            dataPath,
            configDataObj,
            excludeString=excludeString,
            includeString=includeString,
        )
    else:
        print("!! No file numbers found to retrieve sample data")
        pass


def allDataGrabber(configDataObj):
    """Retrieves all data (dark current, background, and sample) defined in
    configDataObj.

    Attempts to retrieve datasets from PANDA_05, BSDIODES, and
    sample name, for a selected scan number. If not all datasets
    can be retrieved, no datasets for the selected scan number
    are saved, and user is notified. Returns four datasets:
    sampleNames, pandaData, bsdiodesData, and gain (gain setting
    of bsdiodes tetramm).
    """
    try:
        bgData = backgroundDataGrabber(configDataObj)
    except:
        print("Could not retreive background data...")

    try:
        sampleData = sampleDataGrabber(configDataObj)
    except:
        print("!! Could not retreive sample data...")

    return bgData, sampleData


def desmearedDataGrabber(scatteringDataObj, configDataObj, sampleName):
    """
    Read in desmeared data for a specified file (sample) name, in a particular
    directory.

    Directory which holds unsmeared USAXS data is stored in the configDatObj,
    and can be set whenever you wish. Desmeared USAXS data must be in order:
    q, I, Ierr. File must be saved as "sampleName_dsm.txt" (sample name taken
    from nexus files).
    """
    print("Retrieving desmeared USAXS data...")
    filePath = configDataObj.desmearedUSAXSDataPath + "/" + sampleName + "_dsm.txt"
    try:
        desmearedData = np.loadtxt(filePath, skiprows=1, delimiter="\t")
        scatteringDataObj.mData.desmearedQ = desmearedData[:, 0]
        scatteringDataObj.mData.desmearedI = desmearedData[:, 1]
        scatteringDataObj.mData.desmearedIerr = desmearedData[:, 2]
    except:
        print("!! Coudn't find necessary file...")
    return scatteringDataObj


def xyeWriter(dirName, scatteringDataObj, configDataObj):
    """
    Write xye data to a text file in the specified folder within the Processing
    directory.

    Removes nan before exporting data.
    """
    folderPath = configDataObj.dataPath + "/processing/" + dirName
    directoryChecker(folderPath)
    scanNosString = "i22"
    for scanNo in scatteringDataObj.scanNos:
        scanNosString += "-"
        scanNosString += str(scanNo)
    fullPath = (
        folderPath
        + "/"
        + scanNosString
        + ".txt"
    )
    fileChecker(fullPath)

    nanMask = np.isnan(scatteringDataObj.mData.I)

    x = scatteringDataObj.mData.Q[~nanMask]
    y = scatteringDataObj.mData.I[~nanMask]
    xe = scatteringDataObj.mData.QSTD[~nanMask]
    ye = scatteringDataObj.mData.ISTD[~nanMask]

    with open(fullPath, "w") as f:
        f.write("q(" + configDataObj.outputUnit + ")\tI\tISTD\tQSTD\n")
        for ii in range(len(x)):
            f.write("{:.4e}\t{:.4e}\t{:.4e}\n".format(x[ii], y[ii], ye[ii], xe[ii]))

    print("Data written to file: " + fullPath)


def nexusWriter(dirName, scatteringDataObj, configDataObj):
    """
    Write nexus data to a text file in the specified folder within the
    Processing directory.

    Removes nan before exporting data. Only exports merged data, tries to
    export merged Q values, and failing that exports merged yaw values.
    """

    folderPath = configDataObj.dataPath + "/processing/" + dirName
    directoryChecker(folderPath)
    scanNosString = "i22"
    for scanNo in scatteringDataObj.scanNos:
        scanNosString += "-"
        scanNosString += str(scanNo)
    fullPath = (
        folderPath
        + "/"
        + scanNosString
        + ".nxs"
    )
    fileChecker(fullPath)

    nanMask = np.isfinite(scatteringDataObj.mData.I)

    try:
        x = scatteringDataObj.mData.Q[nanMask]
        xe = scatteringDataObj.mData.QSTD[nanMask]
        xLabel = "q"
        xUnit = configDataObj.outputUnit
    except:
        x = scatteringDataObj.mData.yaw[nanMask]
        xe = scatteringDataObj.mData.yawSTD[nanMask]
        xLabel = "yaw"
        xUnit = configDataObj.inputUnit

    y = scatteringDataObj.mData.I[nanMask]
    ye = scatteringDataObj.mData.ISTD[nanMask]
    n = scatteringDataObj.mData.binnedN[nanMask]

    with nxopen(fullPath, "w") as f:
        f["sample"] = NXentry()
        f["sample/title"] = scatteringDataObj.sampleName
        f["sample/scans"] = scatteringDataObj.mData.dataScanNos
        f["sample/bgScans"] = scatteringDataObj.mData.bgScanNos
        f["sample/processingSteps"] = scatteringDataObj.processingSteps

        f["processed"] = NXentry()
        f["processed/result"] = NXdata()
        f["processed/result/" + xLabel] = NXfield(x, name=xLabel, unit=xUnit)
        f["processed/result/data"] = NXfield(y, name="I", unit="counts")
        f["processed/result/errors"] = NXfield(ye, name="errors", unit="counts")
        f["processed/result/" + xLabel + "errors"] = NXfield(
            xe, name=xLabel + "Errors", unit=xUnit
        )
        f["processed/result/n"] = NXfield(n, name="binnedPoints")
        f["processed/result"].nxaxes = f["processed/result/" + xLabel]
        f["processed/result"].nxsignal = f["processed/result/data"]
        f["processed/result"].set_default()

    print("Data written to file: " + fullPath)
