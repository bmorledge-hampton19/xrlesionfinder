import os, subprocess, time, re
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
from benbiohelpers.FileSystemHandling.DirectoryHandling import checkDirs, getIsolatedParentDir
from benbiohelpers.CustomErrors import InvalidPathError
from xrlesionfinder.ProjectManagement.UsefulFileSystemFunctions import getDataDirectory
from xrlesionfinder.AlignmentAndFormatting.FindAdapters import findAdapters as findAdaptersFunc
from typing import List


# Write metadata on the parameters for the alignment, for future reference.
def writeMetadata(rawReadsFilePath, adaptorSeqeuncesFilePath, bowtie2IndexBasenamePath, bowtie2StatsFilePath,
                  bowtie2Version = None, customBowtie2Arguments = None):

    metadataFilePath = os.path.join(os.path.dirname(rawReadsFilePath),".metadata")
    with open(metadataFilePath, 'w') as metadataFile:

        if bowtie2Version is None:
            bowtie2Version = subprocess.check_output(("bowtie2","--version"), encoding=("utf-8"))

        metadataFile.write("Path_to_Index:\n" + bowtie2IndexBasenamePath + "\n\n")
        metadataFile.write("Path_to_Adaptor_Sequences:\n" + adaptorSeqeuncesFilePath + "\n\n")
        metadataFile.write("Bowtie2_Version:\n" + bowtie2Version + "\n")
        metadataFile.write("Bowtie2_Stats:\n")
        with open(bowtie2StatsFilePath, 'r') as bowtie2StatsFile:
            atStats = False
            for line in bowtie2StatsFile:
                if not atStats and line.endswith("reads; of these:\n"): atStats = True
                if atStats: metadataFile.write(line)
            if not atStats: raise InvalidPathError("Malformed bowtie2 stats.")
        metadataFile.write('\n')
        os.remove(bowtie2StatsFilePath)
        if customBowtie2Arguments is not None and customBowtie2Arguments:
            metadataFile.write(f"Custom_Bowtie2_arguments: {customBowtie2Arguments}\n\n")


# For each of the given reads files, run the accompyaning bash script to perform the alignment.
def alignXRSeqReads(rawReadsFilePaths: List[str], adapterSequencesFilePath, bowtie2IndexBasenamePath, alignmentBashScriptFilePath, 
                    readCountsOutputFilePath = None, bowtie2BinaryPath = None, threads = 1, customBowtie2Arguments = '',
                    findAdapters = False, pairedEndAlignment = False):

    readCounts = dict()
    scriptStartTime = time.time()
    totalReadsFiles = len(rawReadsFilePaths)
    currentReadFileNum = 0

    # If performing paired end alignment, find pairs for all the given raw reads files.
    if pairedEndAlignment:
        read1FilePaths = list(); read2FilePaths = list()
        for rawReadsFilePath in rawReadsFilePaths:

            # It's possible we have already assigned this file path as a pair of a previous path. Double check!
            if rawReadsFilePath in read1FilePaths or rawReadsFilePath in read2FilePaths: continue

            # Find the pair and assign each pair to their respective list.
            baseName = rawReadsFilePath.rsplit(".fastq",1)[0]
            if baseName.endswith("R1"):
                read1FilePaths.append(rawReadsFilePath)
                pairedBaseName = baseName[:-1] + '2'
                pairedList = read2FilePaths
            elif baseName.endswith("R2"):
                read2FilePaths.append(rawReadsFilePath)
                pairedBaseName = baseName[:-1] + '1'
                pairedList = read1FilePaths
            else: raise InvalidPathError(rawReadsFilePath, "Given path does not end with \"R1\" or \"R2\" "
                                                           "(prior to file extension; e.g. my_reads_R2.fastq.gz is valid.)")

            pairFound = False
            for fileExtension in (".fastq", ".fastq.gz"):
                pairedFilePath = pairedBaseName + fileExtension
                if os.path.exists(pairedFilePath): 
                    pairedList.append(pairedFilePath)
                    print(f"Found fastq file pair with basename: {os.path.basename(pairedBaseName).rsplit('_R',1)[0]}")
                    pairFound = True
                    break

            if not pairFound: raise InvalidPathError(rawReadsFilePath, "No matching pair found. Expected paired files (in same "
                                                                       "directory) ending in \"R1\" and \"R2\", but only found:")

        rawReadsFilePaths = read1FilePaths

    for i, rawReadsFilePath in enumerate(rawReadsFilePaths):

        # Print information about the current file
        currentReadFileNum += 1
        readsFileStartTime = time.time()
        print()
        if pairedEndAlignment: 
            print(f"Processing file pair with basename {os.path.basename(rawReadsFilePath).rsplit('_R1',1)[0]}")
        else: print("Processing file",os.path.basename(rawReadsFilePath))
        print('(',currentReadFileNum,'/',totalReadsFiles,')', sep = '') 

        # If requested find adapters in the current fastq file(s) based on the given list of adapters.
        if findAdapters and not pairedEndAlignment:
            thisAdapterSequencesFilePath = findAdaptersFunc([rawReadsFilePath], adapterSequencesFilePath)[0]
        elif findAdapters and pairedEndAlignment:
            thisAdapterSequencesFilePath = findAdaptersFunc([read1FilePaths[i], read2FilePaths[i]],
                                                            adapterSequencesFilePath, aggregateOutput = True)[0]
        else: thisAdapterSequencesFilePath = adapterSequencesFilePath

        # Make sure the .tmp directory exists and create a path to the bowtie2 stats file.
        tempDir = os.path.join(os.path.dirname(rawReadsFilePath),".tmp")
        checkDirs(tempDir)
        bowtie2StatsFilePath = os.path.join(tempDir, "bowtie2_stats.txt")

        # Run the alignment script.
        if pairedEndAlignment:
            arguments = ["bash", alignmentBashScriptFilePath, "-1", read1FilePaths[i], "-2", read2FilePaths[i],
                         "-a", thisAdapterSequencesFilePath, "-i", bowtie2IndexBasenamePath,
                         "-t", str(threads), "-c", customBowtie2Arguments]
        else:
            arguments = ["bash", alignmentBashScriptFilePath, "-1", rawReadsFilePath, "-a", thisAdapterSequencesFilePath, 
                         "-i", bowtie2IndexBasenamePath, "-t", str(threads), "-c", customBowtie2Arguments]
        if bowtie2BinaryPath is not None: arguments += ["-b", bowtie2BinaryPath]
        subprocess.run(arguments, check = True)

        # If requested, count the number of reads in the original input file(s).
        if readCountsOutputFilePath is not None:
            print("Counting reads in original reads file(s)...")

            zcatProcess = subprocess.Popen(("zcat", rawReadsFilePath), stdout = subprocess.PIPE)
            readCountProcess = subprocess.Popen(("wc", "-l"), stdin = zcatProcess.stdout, stdout = subprocess.PIPE)
            readCount = readCountProcess.communicate()[0].decode("utf8")
            readCounts[os.path.basename(rawReadsFilePath)] = str( round(int(readCount)/4) )

            if pairedEndAlignment:
                zcatProcess = subprocess.Popen(("zcat", read2FilePaths[i]), stdout = subprocess.PIPE)
                readCountProcess = subprocess.Popen(("wc", "-l"), stdin = zcatProcess.stdout, stdout = subprocess.PIPE)
                readCount = readCountProcess.communicate()[0].decode("utf8")
                readCounts[os.path.basename(read2FilePaths[i])] = str( round(int(readCount)/4) )

        # Output information on time elapsed.
        print(f"Time taken to align reads in this file: {time.time() - readsFileStartTime} seconds")
        print(f"Total time spent aligning across all files: {time.time() - scriptStartTime} seconds")

        # Write the metadata.
        writeMetadata(rawReadsFilePath, thisAdapterSequencesFilePath, bowtie2IndexBasenamePath, 
                      bowtie2StatsFilePath, bowtie2BinaryPath, customBowtie2Arguments)

    # Write the read counts if requested.
    if readCountsOutputFilePath is not None:
        with open(readCountsOutputFilePath, 'w') as readCountsOutputFile:
            for rawReadsFileBasename in readCounts:
                readCountsOutputFile.write(rawReadsFileBasename + ": " + readCounts[rawReadsFileBasename] + '\n')


# Removes trimmed reads file paths from a list of fastq reads file paths.
# Returns the filtered list of file paths. Does not alter the original list.
def removeTrimmedAndTmp(unfilteredReadsFilePaths: List[str]):
    return [filePath for filePath in unfilteredReadsFilePaths if 
            getIsolatedParentDir(filePath) != ".tmp" and not "trimmed.fastq" in filePath and
            re.search("trimmed_..\.fastq", filePath) is None]


def parseArgs(args):
    pass
    # TODO: Implement this


def main():

    # Create a simple dialog for selecting the relevant files.
    with TkinterDialog(workingDirectory=getDataDirectory()) as dialog:
        dialog.createMultipleFileSelector("Raw fastq reads:", 0, ".fastq.gz",
                                        ("Fastq Files", (".fastq.gz", ".fastq")), ("fastq Files", ".fastq"),
                                        additionalFileEndings=[".fastq"])
        dialog.createFileSelector("Bowtie2 Index File (Any):", 1, ("Bowtie2 Index File", ".bt2"))

        dialog.createDropdown("Alignment method:", 2, 0, ("Single-end", "Paired-end"))

        with dialog.createDynamicSelector(3, 0) as adaptorSequencesDS:
            adaptorSequencesDS.initDropdownController("Use default adapters for trimming",
                                                      ("Default", "Custom", "Find Adapters", "None"))
            adaptorSequencesSelector = adaptorSequencesDS.initDisplay("Custom", selectionsID = "customAdapterSequences")
            adaptorSequencesSelector.createFileSelector("Custom Adapters Sequences File:", 0, ("Fasta Files", ".fa"))
            adaptorSequencesSelector = adaptorSequencesDS.initDisplay("Find Adapters", selectionsID = "potentialAdapterSequences")
            adaptorSequencesSelector.createFileSelector("Potential Adapter Sequences File:", 0, ("Fasta Files", ".fa"))

        dialog.createDropdown("How many Threads should be used?", 4, 0, ['1', '2', '3', '4', '5', '6', '7', '8'])

        with dialog.createDynamicSelector(5, 0) as additionalOptions:
            additionalOptions.initDropdownController("Additional Options:", ("Hide", "Show"))
            additionalOptionsDialog = additionalOptions.initDisplay("Show")

            with additionalOptionsDialog.createDynamicSelector(0, 0) as bowtie2BinaryDS:
                bowtie2BinaryDS.initCheckboxController("Choose alternative bowtie2 binary")
                bowtie2BinarySelector = bowtie2BinaryDS.initDisplay(True, selectionsID = "bowtieBinary")
                bowtie2BinarySelector.createFileSelector("bowtie2 binary:", 0, ("Any File", "*"))

            with additionalOptionsDialog.createDynamicSelector(1, 0) as readCountsDS:
                readCountsDS.initCheckboxController("Record initial read counts")
                readCountsSelector = readCountsDS.initDisplay(True, selectionsID = "readCounts")
                readCountsSelector.createFileSelector("Save read counts at:", 0, ("Text File", ".txt"), newFile = True)

            with additionalOptionsDialog.createDynamicSelector(2, 0) as customArgsDS:
                customArgsDS.initDropdownController("Custom Bowtie2 Arguments:", ("None", "From File", "Direct Input"))
                customArgsDS.initDisplay("From File", selectionsID = "customArgs").createFileSelector(
                    "Custom Arguments File:", 0, ("Text File", ".txt")
                )
                customArgsDS.initDisplay("Direct Input", selectionsID = "customArgs").createTextField(
                    "Custom Arguments:", 0, 0, defaultText = ''
                )

    # Get the raw reads files, but make sure that no trimmed reads files have tagged along!
    unfilteredRawReadsFilePaths = dialog.selections.getFilePathGroups()[0]
    filteredRawReadsFilePaths = removeTrimmedAndTmp(unfilteredRawReadsFilePaths)

    bowtie2IndexBasenamePath: str = dialog.selections.getIndividualFilePaths()[0]
    bowtie2IndexBasenamePath = bowtie2IndexBasenamePath.rsplit('.', 2)[0]
    if bowtie2IndexBasenamePath.endswith(".rev"): bowtie2IndexBasenamePath = bowtie2IndexBasenamePath.rsplit('.', 1)[0]

    if dialog.selections.getDropdownSelections()[0] == "Paired-end": pairedEndAlignment = True
    else: pairedEndAlignment = False

    findAdapters = False
    if adaptorSequencesDS.getControllerVar() == "Default":
        adaptorSequencesFilePath = os.path.join(os.path.dirname(__file__), "XR-seq_primers.fa")
    elif adaptorSequencesDS.getControllerVar() == "Custom":
        adaptorSequencesFilePath = dialog.selections.getIndividualFilePaths("customAdapterSequences")[0]
    elif adaptorSequencesDS.getControllerVar() == "Find Adapters":
        adaptorSequencesFilePath = dialog.selections.getIndividualFilePaths("potentialAdapterSequences")[0]
        findAdapters = True
    else: adaptorSequencesFilePath = "NONE"

    threads = int(dialog.selections.getDropdownSelections()[1])

    if readCountsDS.getControllerVar():
        readCountsOutputFilePath = dialog.selections.getIndividualFilePaths("readCounts")[0]
    else: readCountsOutputFilePath = None

    if bowtie2BinaryDS.getControllerVar():
        bowtie2BinaryPath = dialog.selections.getIndividualFilePaths("bowtieBinary")[0]
    else: bowtie2BinaryPath = None

    alignmentBashScriptFilePath = os.path.join(os.path.dirname(__file__),"ParseRawReadsToBed.bash")

    if customArgsDS.getControllerVar() == "None":
        customBowtie2Arguments = ''
    elif customArgsDS.getControllerVar() == "From File":
        customBowtie2ArgsFilePath = dialog.selections.getIndividualFilePaths("customArgs")[0]
        with open(customBowtie2ArgsFilePath, 'r') as customBowtie2ArgsFile:
            customBowtie2Arguments = customBowtie2ArgsFile.readline().strip()
    elif customArgsDS.getControllerVar() == "Direct Input":
        customBowtie2Arguments = dialog.selections.getTextEntries("customArgs")[0]

    alignXRSeqReads(filteredRawReadsFilePaths, adaptorSequencesFilePath, bowtie2IndexBasenamePath, alignmentBashScriptFilePath, 
                    readCountsOutputFilePath, bowtie2BinaryPath, threads, customBowtie2Arguments, findAdapters, pairedEndAlignment)


if __name__ == "__main__": main()