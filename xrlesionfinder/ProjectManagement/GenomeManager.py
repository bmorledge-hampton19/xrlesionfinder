# This script manages the known genomes for xrlesionfinder.
import os
from typing import Dict
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
from benbiohelpers.FileSystemHandling.FastaFileIterator import FastaFileIterator
from benbiohelpers.FileSystemHandling.DirectoryHandling import checkDirs, getIsolatedParentDir
from benbiohelpers.FileSystemHandling.BedToFasta import bedToFasta
from benbiohelpers.CustomErrors import checkIfPathExists, InvalidPathError, UserInputError
from xrlesionfinder.ProjectManagement.UsefulFileSystemFunctions import getDataDirectory, getExternalDataDirectory


class GenomeManagerError(Exception):
    "An error class to be raised when a stored genome fasta files cannot be accessed."

class UnrecognizedGenomeError(GenomeManagerError):
    """
    An error class for when a genome name is given that is not recognized
    (i.e. it is not stored in the genome manager's text file).
    """

    def __init__(self, genomeName: str):
        self.genomeName = genomeName

    def __str__(self):
        return f"Could not find a genome with the name \"{self.genomeName}\"."
    
class MissingGenomeFileError(GenomeManagerError):
    "An error class for when the stored file path for a genome no longer exists."
    def __init__(self, genomeName: str, genomeFastaFilePath: str):
        self.genomeName = genomeName
        self.genomeFastaFilePath = genomeFastaFilePath

    def __str__(self):
        return f"Genome Fasta file for {self.genomeName} was not found at the expected location: {self.genomeFastaFilePath}"


def getGenomeListFilePath():
    "Get the path to the file containing the list of known genomes for xrlesionfinder."
    return os.path.join(getExternalDataDirectory,"genomes.txt")


def getGenomes() -> Dict[str]:
    "Return a dictionary of genome fasta file paths with genome names as keys"
    genomes = dict()
    if os.path.exists(getGenomeListFilePath()):
        with open(getGenomeListFilePath(), 'r') as genomeManagerFile:
            for line in genomeManagerFile:
                genomeName,genomeFilePath = line.strip().split(':')
                genomes[genomeName] = genomeFilePath
    return genomes


def getGenomeFastaFilePath(genomeName):
    "Return the path to given genome's fasta file"
    genomes = getGenomes()
    if genomeName not in genomes: raise UnrecognizedGenomeError(genomeName)
    genomeFastaFilePath = genomes[genomeName]
    if os.path.exists(genomeFastaFilePath): return genomeFastaFilePath
    else: raise MissingGenomeFileError(genomeName, genomeFastaFilePath)


def addGenome(genomeFastaFilePath: str, alias = None):
    """
    Add the given genome fasta file to xrlesionfinder's genome manager file.
    If no alias (the colloquial name for the genome) is given, it is generated from the fasta file's name.
    """

    # Make the given path exists and is valid.
    checkIfPathExists(genomeFastaFilePath)
    if not genomeFastaFilePath.endswith(".fa"):
        raise InvalidPathError(genomeFastaFilePath, postPathMessage = "Expected uncompressed fasta file.")
    
    # If no alias was given, derive it from the fasta file name.
    if alias is None: alias = os.path.basename(genomeFastaFilePath).rsplit('.', 1)[0]
    alias = alias.strip() # Also remove leading and trailing whitespace.
    if ':' in alias: raise UserInputError("Invalid character \":\" used in genome alias.")
    if '\n' in alias: raise UserInputError("Invalid character \"\\n\" used in genome alias.")

    # Add the alias and fasta file path to the dictionary of known genomes, overwriting an old entry if necessary.
    print(f"Adding genome fasta file {genomeFastaFilePath} as {alias}.")
    genomes = getGenomes()
    if alias in genomes and genomes[alias] != genomeFastaFilePath:
        print(f"NOTE: This action overwrites an entry which previously pointed to {genomes[alias]}")
    genomes[alias] = genomeFastaFilePath

    # Rewrite the genome manager file with the updated dictionary.
    with open(getGenomeListFilePath(), 'w') as genomeManagerFile:
        for genomeName in sorted(genomes):
            genomeManagerFile.write(f"{genomeName}:{genomes[genomeName]}\n")


def main():

    # Create a simple dialog for selecting the relevant files.
    with TkinterDialog(workingDirectory=getDataDirectory()) as dialog:
        dialog.createFileSelector("GenomeFastaFile", 0, ("Fasta File", ".fa"))
        with dialog.createDynamicSelector(1, 0, 2) as aliasDynSel:
            aliasDynSel.initCheckboxController("Give custom alias")
            aliasDynSel.initDisplay(1, "CustomAlias").createTextField("Custom Alias", 0, 0, defaultText="My_Genome")