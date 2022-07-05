# This script uses samtools and bedtools to convert a sam file to a bed file.
import os, subprocess, time
from typing import List
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
from xrlesionfinder.ProjectManagement.UsefulFileSystemFunctions import getDataDirectory


# For each of the given reads files, run the accompyaning bash script to perform the alignment.
def trimmedFastqToSam(samFilePaths: List[str]):

    for samFilePath in samFilePaths:

        print(f"Working with {os.path.basename(samFilePath)}")

        bamFilePath = samFilePath.replace("sam","bam")
        bedFilePath = samFilePath.rsplit("sam",1)[0] + "bed"

        # convert to bam
        print("Converting sam to bam...")
        if samFilePath.endswith(".gz"):
            zcatProcess = subprocess.Popen(("zcat", samFilePath), stdout=subprocess.PIPE)
            subprocess.check_call(("samtools", "view", "-b", "-o", bamFilePath), stdin=zcatProcess.stdout)
            zcatProcess.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
        else: subprocess.check_call(("samtools", "view", "-b", "-o", bamFilePath, samFilePath))

        # Gzip the results.
        print("Converting bam to bed...")
        with open(bedFilePath, 'w') as bedFile:
            subprocess.check_call(("bedtools", "bamtobed", "-i", bamFilePath), stdout = bedFile)


def main():

    with TkinterDialog(workingDirectory = getDataDirectory()) as dialog:
        dialog.createMultipleFileSelector("Sam Files:", 0, ".sam", ("Sam Files",(".sam.gz",".sam")),
                                          additionalFileEndings = ".sam.gz")

    trimmedFastqToSam(dialog.selections.getFilePathGroups()[0])


if __name__ == "__main__": main()