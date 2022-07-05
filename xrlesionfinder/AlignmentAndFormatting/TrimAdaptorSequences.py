import os, subprocess
from typing import List
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
from xrlesionfinder.ProjectManagement.UsefulFileSystemFunctions import getDataDirectory


# For each of the given reads files, run the accompanying bash script to perform the alignment.
def trimAdaptorSequences(fastqFilePaths: List[str], adaptorSequencesFilePath):

    trimmingBashScriptFilePath = os.path.join(os.path.dirname(__file__),"TrimAdaptorSequences.bash")

    for fastqFilePath in fastqFilePaths:
        if fastqFilePath.endswith("trimmed.fastq.gz"):
            print(f"{os.path.basename(fastqFilePath)} appears to be already trimmed. Skipping.")
            continue
        else: subprocess.run(("bash", trimmingBashScriptFilePath, fastqFilePath, adaptorSequencesFilePath), check = True)


def main():

    with TkinterDialog(workingDirectory = getDataDirectory()) as dialog:
        dialog.createMultipleFileSelector("Raw fastq reads:", 0, ".fastq.gz", 
                                          ("Gzipped fastq Files", ".fastq.gz"))
        dialog.createFileSelector("Adaptor Sequences:", 1, ("Fasta Files", ".fa"))

    trimAdaptorSequences(dialog.selections.getFilePathGroups()[0], dialog.selections.getIndividualFilePaths()[0])


if __name__ == "__main__": main()