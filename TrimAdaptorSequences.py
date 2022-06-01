import os, subprocess
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import getDataDirectory


# For each of the given reads files, run the accompyaning bash script to perform the alignment.
def trimAdaptorSequences(fastqFilePaths, adaptorSequencesFilePath):
    
    trimmingBashScriptFilePath = os.path.join(os.path.dirname(__file__),"TrimAdaptor.bash")

    for fastqFilePath in fastqFilePaths:
        subprocess.run(("bash", trimmingBashScriptFilePath, fastqFilePath, adaptorSequencesFilePath), check = True)


def main():

    with TkinterDialog(workingDirectory=getDataDirectory()) as dialog:
        dialog.createMultipleFileSelector("Raw fastq reads:", 0, ".fastq.gz", 
                                          ("Gzipped fastq Files", ".fastq.gz"))
        dialog.createFileSelector("Adaptor Sequences:", 1, ("Fasta Files", ".fa"))

    trimAdaptorSequences(dialog.selections.getFilePathGroups()[0], dialog.selections.getIndividualFilePaths()[0])


if __name__ == "__main__": main()