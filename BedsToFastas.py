import os
from typing import List
from mutperiodpy.helper_scripts.UsefulBioinformaticsFunctions import bedToFasta
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import getDataDirectory
from mutperiodpy.Tkinter_scripts.TkinterDialog import TkinterDialog

# Given a list of bed files, simply uses the bedToFasta command to convert them to fasta files.
def bedsToFastas(bedFilePaths: List[str], genomeFilePath):
    
    for bedFilePath in bedFilePaths:
        print("Converting",os.path.basename(bedFilePath))
        bedToFasta(bedFilePath, genomeFilePath, bedFilePath.rsplit('.',1)[0]+".fa")


def main():
    # Create a simple dialog for selecting the bed files and corresponding genome file.
    dialog = TkinterDialog(workingDirectory=getDataDirectory())
    dialog.createMultipleFileSelector("Bed Files:", 0, ".bed", ("Bed Files", ".bed"))
    dialog.createFileSelector("Genome Fasta File:", 1, ("Fasta File", ".fa"))

    dialog.mainloop()
    if dialog.selections is None: quit()

    bedsToFastas(dialog.selections.getFilePathGroups()[0], dialog.selections.getIndividualFilePaths()[0])


if __name__ == "__main__": main()