import subprocess
from mutperiodpy.Tkinter_scripts.TkinterDialog import TkinterDialog
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import getDataDirectory
from typing import List
import os, subprocess

# For each given file, create a new file where every read position is unique.  Reads with duplicate locations are trimmed down to the one read.
def trimDuplicateReads(readsFilePaths: List[str]):

    for readsFilePath in readsFilePaths:

        print("\nWorking in",os.path.basename(readsFilePath))

        # Create a file to output results to that doesn't include duplicates
        noDupsFilePath = readsFilePath.rsplit('.',1)[0] + "_no_dups.bed"

        # Sort the file first, just to be sure!
        print("Sorting...")
        subprocess.run( ("sort", "-k1,1", "-k2,3n", "-o", readsFilePath, readsFilePath), check = True)

        encounteredReads = dict()

        # Iterate through the sorted reads file, writing each line to the new file but omitting any duplicate reads beyond the first.
        print("Removing excess duplicates...")
        with open(readsFilePath, 'r') as readsFile:
            with open(noDupsFilePath, 'w') as noDupsFile:

                lastReadLocation = ''

                for line in readsFile:

                    readLocation = ' '.join(line.split()[:3])

                    # Write this read only if it does not match the previous read.
                    if readLocation != lastReadLocation:

                        noDupsFile.write(line)

                        assert readLocation not in encounteredReads, "Read has been encountered before: " + readLocation
                        encounteredReads[readLocation] = None

                        lastReadLocation = readLocation


def main():

    # Create a simple dialog for selecting the gene designation files.
    dialog = TkinterDialog(workingDirectory=getDataDirectory())
    dialog.createMultipleFileSelector("Bed formatted reads:", 0, ".bed", 
                                      ("Bed Files", ".bed"))

    dialog.mainloop()

    if dialog.selections is None: quit()

    trimDuplicateReads(dialog.selections.getFilePathGroups()[0])

if __name__ == "__main__": main()