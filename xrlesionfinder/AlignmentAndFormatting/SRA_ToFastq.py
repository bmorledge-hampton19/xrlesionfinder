# This script takes a list of SRA run accession IDs and corresponding names 
# (both formatted as simple newline-separated text files),
# and uses sra-tools to download and convert the reads to fastq.gz format.

import os, subprocess, time, shutil
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
from benbiohelpers.FileSystemHandling.DirectoryHandling import checkDirs
from xrlesionfinder.ProjectManagement.UsefulFileSystemFunctions import getDataDirectory

# Given a file path to a list of run accession IDs and an optional file path to their corresponding names,
# uses sra-tools to retrieve the fastq sequences.
def sRA_ToFastq(runAccessionIDsFilePath, namesFilePath = None):
    
    # Get the run accession IDs from the input file
    with open(runAccessionIDsFilePath, 'r') as runAccessionIDsFile:
        runAccessionIDs = [line.strip() for line in runAccessionIDsFile]

    # Determine the names for the output files.
    if namesFilePath is None: names = runAccessionIDs
    else: 
        with open(namesFilePath, 'r') as namesFile:
            names = [line.strip() for line in namesFile]

    startTime = time.time()
    gzippedFastqFiles = list()
    for runAccessionID, name in zip(runAccessionIDs, names):

        print(f"\nRetrieving reads for {runAccessionID}")
        

        # Generate the paths for output.
        alignmentFilesDir = os.path.join(os.path.dirname(runAccessionIDsFilePath), name, "alignment_files")
        checkDirs(alignmentFilesDir)
        fastqOutputFilePath = os.path.join(alignmentFilesDir, name + ".fastq")

        # Run prefetch.
        prefetchStartTime = time.time()
        print("\nRunning prefetch...")
        subprocess.check_call(("prefetch", "-p", "-O", alignmentFilesDir, runAccessionID))
        print(f"Time to fetch: {time.time() - prefetchStartTime} seconds")

        # Run fasterqdump.
        fasterqDumpStartTime = time.time()
        print("\nRunning fasterq-dump...")
        subprocess.check_call(("fasterq-dump", "-p", "-o", fastqOutputFilePath, runAccessionID), cwd = alignmentFilesDir)
        print(f"Time to retrieve fastq: {time.time() - fasterqDumpStartTime} seconds")

        # gzip the results.
        print("\ngzipping results...")
        gzipStartTime = time.time()
        for item in os.listdir(alignmentFilesDir):
            if item.endswith(".fastq"): 
                print(f"gzipping {item}")
                subprocess.check_call(("gzip", "-f", os.path.join(alignmentFilesDir, item)))
                gzippedFastqFiles.append(os.path.join(alignmentFilesDir, item+".gz"))
        print(f"Time to gzip files: {time.time() - gzipStartTime} seconds")

        # Delete prefetched directory
        print("\nDeleting prefetched directory...")
        shutil.rmtree(os.path.join(alignmentFilesDir,runAccessionID))

        print(f"**Total time so far: {time.time() - startTime} seconds**")

    return gzippedFastqFiles


def main():

    # Create the UI.
    with TkinterDialog(workingDirectory = getDataDirectory()) as dialog:
        dialog.createFileSelector("SRA run accession IDs:", 0, ("text file",".txt"))
        with dialog.createDynamicSelector(1, 0) as namesDynSel:
            namesDynSel.initCheckboxController("Specify Custom Names:")
            namesDynSel.initDisplay(True, "namesFilePath").createFileSelector(
                "Names File:", 0, ("text file",".txt")
            )

    # Retrieve values from user input.
    selections = dialog.selections
    if namesDynSel.getControllerVar(): namesFilePath = selections.getIndividualFilePaths("namesFilePath")[0]
    else: namesFilePath = None

    sRA_ToFastq(selections.getIndividualFilePaths()[0], namesFilePath)
    

if __name__ == "__main__": main()