import os, subprocess
from mutperiodpy.Tkinter_scripts.TkinterDialog import TkinterDialog
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import getDataDirectory


# Write metadata on the parameters for the alignment, for future reference.
def writeMetadata(rawReadsFilePath, adaptorSeqeuncesFilePath, bowtie2IndexBasenamePath):

    metadataFilePath = os.path.join(os.path.dirname(rawReadsFilePath),".metadata")
    with open(metadataFilePath, 'w') as metadataFile:

        bowtie2Version = subprocess.check_output(("bowtie2","--version"), encoding=("utf-8"))

        metadataFile.write("Path_to_Index:\n" + bowtie2IndexBasenamePath + "\n\n")
        metadataFile.write("Path_to_Adaptor_Sequences:\n" + adaptorSeqeuncesFilePath + "\n\n")
        metadataFile.write("Bowtie2_Version: " + bowtie2Version + "\n\n")


# For each of the given reads files, run the accompyaning bash script to perform the alignment.
def alignXRSeqReads(rawReadsFilePaths, adaptorSequencesFilePath, bowtie2IndexBasenamePath, alignmentBashScriptFilePath, 
                    readCountsOutputFilePath = None, bowtie2BinaryPath = None):
    
    readCounts = dict()
    totalReadsFiles = len(rawReadsFilePaths)
    currentReadFileNum = 0
    for rawReadsFilePath in rawReadsFilePaths:

        # Print information about the current file
        currentReadFileNum += 1
        print()
        print("Processing file",os.path.basename(rawReadsFilePath))
        print('(',currentReadFileNum,'/',totalReadsFiles,')', sep = '') 

        # Run the alignment script.
        arguments = ["bash", alignmentBashScriptFilePath, rawReadsFilePath, adaptorSequencesFilePath, bowtie2IndexBasenamePath]
        if bowtie2BinaryPath is not None: arguments.append(bowtie2BinaryPath)
        subprocess.run(arguments, check = True)

        # If requested, count the number of reads in the original input file.
        if readCountsOutputFilePath is not None:
            print("Counting reads in original reads file...")

            zcatProcess = subprocess.Popen(("zcat", rawReadsFilePath), stdout = subprocess.PIPE)
            readCountProcess = subprocess.Popen(("wc", "-l"), stdin = zcatProcess.stdout, stdout = subprocess.PIPE)
            readCount = readCountProcess.communicate()[0].decode("utf8")
            readCounts[os.path.basename(rawReadsFilePath)] = str( (int(readCount)-1)/4 )

        # Write the metadata.
        writeMetadata(rawReadsFilePath, adaptorSequencesFilePath, bowtie2IndexBasenamePath)

    # Write the read counts if requested.
    if readCountsOutputFilePath is not None:
        with open(readCountsOutputFilePath, 'w') as readCountsOutputFile:
            for rawReadsFileBasename in readCounts:
                readCountsOutputFile.write(rawReadsFileBasename + ": " + readCounts[rawReadsFileBasename] + '\n')

def main():

    # Create a simple dialog for selecting the gene designation files.
    dialog = TkinterDialog(workingDirectory=getDataDirectory())
    dialog.createMultipleFileSelector("Raw fastq reads:", 0, ".fastq.gz", 
                                      ("Gzipped fastq Files", ".fastq.gz"))
    dialog.createFileSelector("Bowtie2 Index File (Any):", 1, ("Bowtie2 Index File", ".bt2"))

    adaptorSequencesDS = dialog.createDynamicSelector(2, 0)
    adaptorSequencesDS.initCheckboxController("Specify adaptor sequences to trim")
    adaptorSequencesSelector = adaptorSequencesDS.initDisplay(True, selectionsID = "adaptorSequences")
    adaptorSequencesSelector.createFileSelector("Adaptor Sequences:", 0, ("Fasta Files", ".fa"))
    adaptorSequencesDS.initDisplayState()

    bowtie2BinaryDS = dialog.createDynamicSelector(3, 0)
    bowtie2BinaryDS.initCheckboxController("Choose alternative bowtie2 binary")
    bowtie2BinarySelector = bowtie2BinaryDS.initDisplay(True, selectionsID = "bowtieBinary")
    bowtie2BinarySelector.createFileSelector("bowtie2 binary:", 0, ("Any File", "*"))
    bowtie2BinaryDS.initDisplayState()
    
    readCountsDS = dialog.createDynamicSelector(4, 0)
    readCountsDS.initCheckboxController("Record initial read counts")
    readCountsSelector = readCountsDS.initDisplay(True, selectionsID = "readCounts")
    readCountsSelector.createFileSelector("Save read counts at:", 0, ("Text File", ".txt"), newFile = True)
    readCountsDS.initDisplayState()

    dialog.mainloop()

    if dialog.selections is None: quit()

    # Get the raw reads files, but make sure that no trimmed reads files have tagged along!
    unfilteredRawReadsFilePaths = dialog.selections.getFilePathGroups()[0]
    filteredRawReadsFilePaths = list()
    for unfilteredRawReadsFilePath in unfilteredRawReadsFilePaths:
        if not unfilteredRawReadsFilePath.endswith("trimmed.fastq.gz"):
            filteredRawReadsFilePaths.append(unfilteredRawReadsFilePath)

    bowtie2IndexBasenamePath: str = dialog.selections.getIndividualFilePaths()[0]
    bowtie2IndexBasenamePath = bowtie2IndexBasenamePath.rsplit('.', 2)[0]
    if bowtie2IndexBasenamePath.endswith(".rev"): bowtie2IndexBasenamePath = bowtie2IndexBasenamePath.rsplit('.', 1)[0]

    if adaptorSequencesDS.getControllerVar():
        adaptorSequencesFilePath = dialog.selections.getIndividualFilePaths("adaptorSequences")[0]
    else: adaptorSequencesFilePath = "NONE"

    if readCountsDS.getControllerVar():
        readCountsOutputFilePath = dialog.selections.getIndividualFilePaths("readCounts")[0]
    else: readCountsOutputFilePath = None

    if bowtie2BinaryDS.getControllerVar():
        bowtie2BinaryPath = dialog.selections.getIndividualFilePaths("bowtieBinary")[0]
    else: bowtie2BinaryPath = None

    alignmentBashScriptFilePath = os.path.join(os.path.dirname(__file__),"ParseRawReadsToBed.bash")

    alignXRSeqReads(filteredRawReadsFilePaths, adaptorSequencesFilePath, bowtie2IndexBasenamePath, alignmentBashScriptFilePath, 
                    readCountsOutputFilePath, bowtie2BinaryPath)


if __name__ == "__main__": main()