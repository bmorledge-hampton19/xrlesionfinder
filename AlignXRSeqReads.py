import os, subprocess
from mutperiodpy.Tkinter_scripts.TkinterDialog import TkinterDialog
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import getDataDirectory

# For each of the given reads files, run the accompyaning bash script to perform the alignment.
def alignXRSeqReads(rawReadsFilePaths, adapterSequencesFilePath, bowtie2IndexBasenamePath, alignmentBashScriptFilePath, readCountsOutputFilePath = None):
    
    readCounts = dict()
    totalReadsFiles = len(rawReadsFilePaths)
    currentReadFileNum = 0
    for rawReadsFilePath in rawReadsFilePaths:

        # Print information about the current file
        currentReadFileNum += 1
        print()
        print("Processing file",os.path.basename(rawReadsFilePath))
        print("(",currentReadFileNum,'/',totalReadsFiles,')', sep = '') 

        # Run the alignment script.
        subprocess.run(("bash", alignmentBashScriptFilePath, rawReadsFilePath, adapterSequencesFilePath, bowtie2IndexBasenamePath), 
                       check = True)

        # If requested, count the number of reads in the original input file.
        if readCountsOutputFilePath is not None:
            print("Counting reads in original reads file...")

            zcatProcess = subprocess.Popen(("zcat", rawReadsFilePath), stdout = subprocess.PIPE)
            readCountProcess = subprocess.Popen(("wc", "-l"), stdin = zcatProcess.stdout, stdout = subprocess.PIPE)
            readCount = readCountProcess.communicate()[0].decode("utf8")
            readCounts[os.path.basename(rawReadsFilePath)] = str( (int(readCount)-1)/4 )

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
    dialog.createFileSelector("Adapter Sequences:", 1, ("Fasta Files", ".fa"))
    dialog.createFileSelector("Bowtie2 Index File (Any):", 2, ("Bowtie2 Index File", ".bt2"))
    readCountsDS = dialog.createDynamicSelector(3, 0)
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

    adapterSequencesFilePath = dialog.selections.getIndividualFilePaths()[0]

    bowtie2IndexBasenamePath: str = dialog.selections.getIndividualFilePaths()[1]
    bowtie2IndexBasenamePath = bowtie2IndexBasenamePath.rsplit('.', 2)[0]
    if bowtie2IndexBasenamePath.endswith(".rev"): bowtie2IndexBasenamePath = bowtie2IndexBasenamePath.rsplit('.', 1)[0]

    if readCountsDS.getControllerVar():
        readCountsOutputFilePath = dialog.selections.getIndividualFilePaths("readCounts")[0]
    else: readCountsOutputFilePath = None


    alignmentBashScriptFilePath = os.path.join(os.path.dirname(__file__),"ParseRawReadsToBed.bash")

    alignXRSeqReads(filteredRawReadsFilePaths, adapterSequencesFilePath, bowtie2IndexBasenamePath, alignmentBashScriptFilePath, readCountsOutputFilePath)


if __name__ == "__main__": main()