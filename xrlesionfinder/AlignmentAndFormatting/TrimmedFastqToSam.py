# This script isolates the alignment process from the AlignXRSeqReads script.
import os, subprocess, time
from typing import List
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
from xrlesionfinder.ProjectManagement.UsefulFileSystemFunctions import getDataDirectory


# For each of the given reads files, run the accompyaning bash script to perform the alignment.
def trimmedFastqToSam(fastqFilePaths: List[str], bowtie2IndexBasenamePath, bowtie2BinaryPath, 
                      threads, customBowtie2Arguments):

    scriptStartTime = time.time()
    totalReadsFiles = len(fastqFilePaths)
    currentReadFileNum = 0

    for fastqFilePath in fastqFilePaths:

        # Get things set up...
        currentReadFileNum += 1
        readsFileStartTime = time.time()
        print(f"\nWorking with {os.path.basename(fastqFilePath)}")
        print('(',currentReadFileNum,'/',totalReadsFiles,')', sep = '') 
        if bowtie2BinaryPath is None: bowtie2BinaryPath = "bowtie2"
        samOutputFilePath = fastqFilePath.rsplit("_trimmed.fastq", 1)[0] + ".sam"

        # Perform the alignment
        print("Aligning reads using bowtie2...")
        subprocess.run((bowtie2BinaryPath, "-x", bowtie2IndexBasenamePath, "-U", fastqFilePath,
                        "-S", samOutputFilePath, "-p", threads, customBowtie2Arguments), check = True)

        # Gzip the results.
        print("gzipping results...")
        subprocess.run(("gzip", "-f", samOutputFilePath))

        # Output information on time elapsed.
        print(f"Time taken to process this file: {time.time() - readsFileStartTime} seconds")
        print(f"Total time spent processing across all files: {time.time() - scriptStartTime} seconds")


def main():

    with TkinterDialog(workingDirectory = getDataDirectory()) as dialog:
        dialog.createMultipleFileSelector("Trimmed fastq reads:", 0, "trimmed.fastq.gz", 
                                          ("Gzipped fastq Files", ".fastq.gz"))
        dialog.createFileSelector("Bowtie2 Index File (Any):", 1, ("Bowtie2 Index File", ".bt2"))

        with dialog.createDynamicSelector(2, 0) as bowtie2BinaryDS:
            bowtie2BinaryDS.initCheckboxController("Choose alternative bowtie2 binary")
            bowtie2BinarySelector = bowtie2BinaryDS.initDisplay(True, selectionsID = "bowtieBinary")
            bowtie2BinarySelector.createFileSelector("bowtie2 binary:", 0, ("Any File", "*"))

        dialog.createDropdown("Number of CPU threads to use:", 3, 0, ('1', '2', '3', '4'))

        with dialog.createDynamicSelector(4, 0) as customArgsDS:
            customArgsDS.initDropdownController("Custom Bowtie2 Arguments:", ("None", "From File", "Direct Input"))
            customArgsDS.initDisplay("From File", selectionsID = "customArgs").createFileSelector(
                "Custom Arguments File:", 0, ("Text File", ".txt")
            )
            customArgsDS.initDisplay("Direct Input", selectionsID = "customArgs").createTextField(
                "Custom Arguments:", 0, 0, defaultText = ''
            )

    fastqFilePaths = dialog.selections.getFilePathGroups()[0]

    bowtie2IndexBasenamePath: str = dialog.selections.getIndividualFilePaths()[0]
    bowtie2IndexBasenamePath = bowtie2IndexBasenamePath.rsplit('.', 2)[0]
    if bowtie2IndexBasenamePath.endswith(".rev"): bowtie2IndexBasenamePath = bowtie2IndexBasenamePath.rsplit('.', 1)[0]

    if bowtie2BinaryDS.getControllerVar():
        bowtie2BinaryPath = dialog.selections.getIndividualFilePaths("bowtieBinary")[0]
    else: bowtie2BinaryPath = None

    threads = dialog.selections.getDropdownSelections()[0]
    
    if customArgsDS.getControllerVar() == "None":
        customBowtie2Arguments = ''
    elif customArgsDS.getControllerVar() == "From File":
        customBowtie2ArgsFilePath = dialog.selections.getIndividualFilePaths("customArgs")[0]
        with open(customBowtie2ArgsFilePath, 'r') as customBowtie2ArgsFile:
            customBowtie2Arguments = customBowtie2ArgsFile.readline().strip()
    elif customArgsDS.getControllerVar() == "Direct Input":
        customBowtie2Arguments = dialog.selections.getTextEntries("customArgs")[0]

    trimmedFastqToSam(fastqFilePaths, bowtie2IndexBasenamePath,
                      bowtie2BinaryPath, threads, customBowtie2Arguments)


if __name__ == "__main__": main()