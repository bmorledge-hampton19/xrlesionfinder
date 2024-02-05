# This script serves as a quick "all-in-one" pipeline to call lesion positions from XR-seq data.
# Inputs include the XR-seq data (in sam or fastq format; a basic XR-seq alignment protocol will be invoked if fastq reads are given),
# a genome version, a list of valid read lengths, and information about the lesion in question.
import os
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
from xrlesionfinder.ProjectManagement.UsefulFileSystemFunctions import getDataDirectory
from xrlesionfinder.ProjectManagement.GenomeManager import getGenomes


def fullLesionFinder(inputFilePaths, associatedGenome, validReadLengths, lesionMismatchSignature):
    pass

    # If given 





def main():

    with TkinterDialog(workingDirectory = getDataDirectory(), title = "Quick XR-seq Lesion Finder") as dialog:
        dialog.createMultipleFileSelector("Input files:", 0, ".sam", ("Sam Files",(".sam", ".sam.gz")), ("Fastq Files", (".fastq", ".fastq.gz")),
                                          additionalFileEndings = ".sam.gz")
        dialog.createDropdown("Genome:", 1, 0, getGenomes().keys())
        with dialog.createDynamicSelector(2, 0, 2) as lesionInfoDynSel:
            lesionInfoDynSel.initDropdownController("Lesion Mismatch Information:", ("Use Presets", "Give Custom Signature"))



if __name__ == "__main__": main()