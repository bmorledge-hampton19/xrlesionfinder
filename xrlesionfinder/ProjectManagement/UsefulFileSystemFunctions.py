# This script contains various functions that I think will often be useful when managing filesystems for this project.

import os
from benbiohelpers.FileSystemHandling.DirectoryHandling import checkDirs
from benbiohelpers.CustomErrors import UserInputError, InvalidPathError

# Get the data directory for analyzeXR-seq, creating it from user input if necessary.
def getDataDirectory():

    # Check for the text file which should contain the path to the data directory.
    dataDirectoryTextFilePath = os.path.join(os.getenv("HOME"), ".xrlesionfinder", "data_dir.txt")

    # If it exists, return the directory path within.
    if os.path.exists(dataDirectoryTextFilePath):
        with open(dataDirectoryTextFilePath, 'r') as dataDirectoryTextFile:
            
            dataDirectory = dataDirectoryTextFile.readline().strip()
            
            # Double check to make sure the data directory is still intact.  
            # If it isn't, inform the user, and progress through the function to recreate it.
            if not os.path.isdir(dataDirectory):
                print("Data directory not found at expected location: {}".format(dataDirectory))
                print("Please select a new location to create a data directory.")
            else: return dataDirectory

    # Create a simple dialog to select a new data directory location.
    # NOTE: The following code is not part of an else statement because the above "if" block will return
    # the data directory if it proceeds correctly, and if it doesn't, the data directory text file
    # needs to be recreated anyway.
    from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog, Selections
    checkDirs(os.path.dirname(dataDirectoryTextFilePath))
    dialog = TkinterDialog(workingDirectory = os.path.dirname(dataDirectoryTextFilePath))
    dialog.createFileSelector("Location to create new data directory:",0,("Fasta Files",".fa"), directory = True)

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    selections: Selections = dialog.selections
    dataDirectoryDirectory = selections.getIndividualFilePaths()[0]

    # Make sure a valid, writeable directory was given.  Then create the new directory (if it doesn't exist already), 
    # write it to the text file, and return it!  (Also create the __external_data directory.)
    if not os.path.exists(dataDirectoryDirectory): 
        raise UserInputError("Given directory: " + dataDirectoryDirectory + " does not exist.")

    dataDirectory = os.path.join(dataDirectoryDirectory,"xrlesionfinder_data")
    try:
        checkDirs(dataDirectory)
        checkDirs(os.path.join(dataDirectory,"__external_data"))
    except IOError:
        raise InvalidPathError(dataDirectoryDirectory, "Given location for data directory is not writeable:")
    with open(dataDirectoryTextFilePath, 'w') as dataDirectoryTextFile:
        dataDirectoryTextFile.write(dataDirectory + '\n')
    getExternalDataDirectory()
    return dataDirectory


# Get the external data directory, creating it if necessary.
def getExternalDataDirectory(): 
    
    externalDataDirectory = os.path.join(getDataDirectory(), "__external_data")
    checkDirs(externalDataDirectory)
    return externalDataDirectory