import os, subprocess, time
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
from benbiohelpers.FileSystemHandling.DirectoryHandling import checkDirs, getIsolatedParentDir
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
            metadataFile.write(''.join(bowtie2StatsFile.readlines()))
        os.remove(bowtie2StatsFilePath)
        if customBowtie2Arguments is not None and customBowtie2Arguments:
            metadataFile.write(f"Custom_Bowtie2_arguments: {customBowtie2Arguments}\n\n")


# For each of the given reads files, run the accompyaning bash script to perform the alignment.
def alignXRSeqReads(rawReadsFilePaths, adapterSequencesFilePath, bowtie2IndexBasenamePath, alignmentBashScriptFilePath, 
                    readCountsOutputFilePath = None, bowtie2BinaryPath = None, threads = 1, customBowtie2Arguments = '',
                    findAdapters = False):

    readCounts = dict()
    scriptStartTime = time.time()
    totalReadsFiles = len(rawReadsFilePaths)
    currentReadFileNum = 0
    for rawReadsFilePath in rawReadsFilePaths:

        # Print information about the current file
        currentReadFileNum += 1
        readsFileStartTime = time.time()
        print()
        print("Processing file",os.path.basename(rawReadsFilePath))
        print('(',currentReadFileNum,'/',totalReadsFiles,')', sep = '') 

        # If requested find adapters in the current fastq file based on the given list of adapters.
        if findAdapters: thisAdapterSequencesFilePath = findAdaptersFunc([rawReadsFilePath], adapterSequencesFilePath)[0]
        else: thisAdapterSequencesFilePath = adapterSequencesFilePath

        # Make sure the .tmp directory exists and create a path to the bowtie2 stats file.
        tempDir = os.path.join(os.path.dirname(rawReadsFilePath),".tmp")
        checkDirs(tempDir)
        bowtie2StatsFilePath = os.path.join(tempDir, "bowtie2_stats.txt")

        # Run the alignment script.
        arguments = ["bash", alignmentBashScriptFilePath, rawReadsFilePath, thisAdapterSequencesFilePath, 
                     bowtie2IndexBasenamePath, str(threads), customBowtie2Arguments]
        if bowtie2BinaryPath is not None: arguments.append(bowtie2BinaryPath)
        subprocess.run(arguments, check = True)

        # If requested, count the number of reads in the original input file.
        if readCountsOutputFilePath is not None:
            print("Counting reads in original reads file...")

            zcatProcess = subprocess.Popen(("zcat", rawReadsFilePath), stdout = subprocess.PIPE)
            readCountProcess = subprocess.Popen(("wc", "-l"), stdin = zcatProcess.stdout, stdout = subprocess.PIPE)
            readCount = readCountProcess.communicate()[0].decode("utf8")
            readCounts[os.path.basename(rawReadsFilePath)] = str( (int(readCount)-1)/4 )

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
            getIsolatedParentDir(filePath) != ".tmp" and not "trimmed.fastq" in filePath]


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

        with dialog.createDynamicSelector(2, 0) as adaptorSequencesDS:
            adaptorSequencesDS.initDropdownController("Use default adapters for trimming", ("Default", "Custom", "Find Adapters", "None"))
            adaptorSequencesSelector = adaptorSequencesDS.initDisplay("Custom", selectionsID = "customAdapterSequences")
            adaptorSequencesSelector.createFileSelector("Custom Adapters Sequences File:", 0, ("Fasta Files", ".fa"))
            adaptorSequencesSelector = adaptorSequencesDS.initDisplay("Find Adapters", selectionsID = "potentialAdapterSequences")
            adaptorSequencesSelector.createFileSelector("Potential Adapter Sequences File:", 0, ("Fasta Files", ".fa"))

        dialog.createDropdown("How many Threads should be used?", 3, 0, ['1', '2', '3', '4', '5', '6', '7', '8'])

        with dialog.createDynamicSelector(4, 0) as bowtie2BinaryDS:
            bowtie2BinaryDS.initCheckboxController("Choose alternative bowtie2 binary")
            bowtie2BinarySelector = bowtie2BinaryDS.initDisplay(True, selectionsID = "bowtieBinary")
            bowtie2BinarySelector.createFileSelector("bowtie2 binary:", 0, ("Any File", "*"))

        with dialog.createDynamicSelector(5, 0) as readCountsDS:
            readCountsDS.initCheckboxController("Record initial read counts")
            readCountsSelector = readCountsDS.initDisplay(True, selectionsID = "readCounts")
            readCountsSelector.createFileSelector("Save read counts at:", 0, ("Text File", ".txt"), newFile = True)

        with dialog.createDynamicSelector(6, 0) as customArgsDS:
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

    findAdapters = False
    if adaptorSequencesDS.getControllerVar() == "Default":
        adaptorSequencesFilePath = os.path.join(os.path.dirname(__file__), "XR-seq_primers.fa")
    elif adaptorSequencesDS.getControllerVar() == "Custom":
        adaptorSequencesFilePath = dialog.selections.getIndividualFilePaths("customAdapterSequences")[0]
    elif adaptorSequencesDS.getControllerVar() == "Find Adapters":
        adaptorSequencesFilePath = dialog.selections.getIndividualFilePaths("potentialAdapterSequences")[0]
        findAdapters = True
    else: adaptorSequencesFilePath = "NONE"

    threads = int(dialog.selections.getDropdownSelections()[0])

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
                    readCountsOutputFilePath, bowtie2BinaryPath, threads, customBowtie2Arguments, findAdapters)


if __name__ == "__main__": main()