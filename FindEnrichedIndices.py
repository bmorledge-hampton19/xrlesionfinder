import os
from enum import Enum
from typing import Dict, List
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import getDataDirectory
from benbiohelpers.FileSystemHandling.FastaFileIterator import FastaFileIterator
from benbiohelpers.FileSystemHandling.DirectoryHandling import checkDirs, getIsolatedParentDir
from benbiohelpers.FileSystemHandling.BedToFasta import bedToFasta


class BaseFrequencyTable:

    class TrackedFeature(Enum):

        singleBase = 1
        dipys = 2

    def __init__(self, fromStartValues, fromEndValues, trackedFeature):

        # Store given features
        self.fromStartValues = fromStartValues
        self.fromEndValues = fromEndValues
        assert trackedFeature in self.TrackedFeature, "Unrecognized Tracked Feature: " + str(trackedFeature)
        self.trackedFeature = trackedFeature

    
    # Generates the dictionary of base frequencies.
    def generateBaseFrequencyTable(self, sequences):

        # Create the dictionary which will contain the frequencies for each base, at each given location.
        # The first dictionary uses the "tracked feature" as its key.
        # The second dictionary uses int values as a key.  Positive values represent "fromStart" values, 
        # and negative values represent "fromEnd values."  Both are 1-based.
        self.baseFrequencies: Dict[str, Dict[int, float]] = dict()

        if self.trackedFeature == self.TrackedFeature.singleBase:
            for base in ('A','C','G','T'):
                self.baseFrequencies[base] = dict()
        elif self.trackedFeature == self.TrackedFeature.dipys:
            for dipy in ("CC","CT","TT","TC","dipys"):
                self.baseFrequencies[dipy] = dict()

        for feature in self.baseFrequencies:
            for fromStartValue in self.fromStartValues: self.baseFrequencies[feature][fromStartValue] = 0
            for fromEndValue in self.fromEndValues: self.baseFrequencies[feature][-fromEndValue] = 0

        # Iterate through all the sequences, counting features at each .
        for sequence in sequences:

            # Record all features at requested start values.
            for fromStartValue in self.fromStartValues:

                if self.trackedFeature == self.TrackedFeature.singleBase:
                    feature = sequence[fromStartValue - 1]
                elif self.trackedFeature == self.TrackedFeature.dipys:
                    feature = sequence[fromStartValue - 1:fromStartValue + 1]
                
                # This check might be important if we encounter "N"
                if feature in self.baseFrequencies: self.baseFrequencies[feature][fromStartValue] += 1

            for fromEndValue in self.fromEndValues:

                if self.trackedFeature == self.TrackedFeature.singleBase:
                    feature = sequence[-fromEndValue]
                elif self.trackedFeature == self.TrackedFeature.dipys:
                    # Edge case here because [-2:0] returns nothing, unlike [-2:]
                    if -fromEndValue == -1: feature = sequence[-fromEndValue-1:]
                    else: feature = sequence[-fromEndValue - 1:-fromEndValue + 1]
                
                # Since we're only counting dipys, make sure this is one of those before counting it!
                if feature in self.baseFrequencies: self.baseFrequencies[feature][-fromEndValue] += 1

        # Convert the raw counts to frequencies.
        sequenceNum = len(sequences)
        if sequenceNum > 0:
            for feature in self.baseFrequencies:
                for pos in self.baseFrequencies[feature]:
                    self.baseFrequencies[feature][pos] = self.baseFrequencies[feature][pos] / sequenceNum

        # If the tracked feature is dipys, add a new feature that is the sum of the 4 dipys.
        if self.trackedFeature == self.TrackedFeature.dipys:
            for pos in self.baseFrequencies["CC"]:
                self.baseFrequencies["dipys"][pos] = sum(self.baseFrequencies[dipy][pos] for dipy in ("CC","CT","TT","TC"))
    

    # Given a feature, return a tuple of the frequency of that feature and the position it is present at.
    def getMaxFrequencyAndPos(self, feature, getSecondPlace = False):

        if feature == self.TrackedFeature.dipys: feature = "dipys"

        assert feature in self.baseFrequencies, "Requested feature, \"" + feature + "\", is not available for this base frequency table."

        maxFrequencyPos = max(self.baseFrequencies[feature], key = self.baseFrequencies[feature].get)

        # If requested, find second place instead.
        if getSecondPlace:
            sansMaxBaseFrequencies = self.baseFrequencies[feature].copy()
            sansMaxBaseFrequencies.pop(maxFrequencyPos)
            maxFrequencyPos = max(sansMaxBaseFrequencies, key = self.baseFrequencies[feature].get)

        maxFrequency = self.baseFrequencies[feature][maxFrequencyPos]
        if maxFrequency == 0: maxFrequencyPos = 0
        return (maxFrequency, self.formatPos(maxFrequencyPos))


    # For tracked dipys, converts the position to a "half-base" position.  For single base positions, just returns the original position.
    def formatPos(self, dipyPos):
        
        if dipyPos == 0: return None

        if self.trackedFeature == self.TrackedFeature.dipys:
            if dipyPos < 0: return dipyPos - 0.5
            else: return dipyPos + 0.5
        else: return dipyPos


    def getBaseFrequencies(self): return self.baseFrequencies


# Reads the given fasta file line by line and filters out any reads of inappropriate length.  
# Returns a dictionary of lists of sequences by length
def getReadSequencesByLength(fastaFilePath, readSizeRange):

    # Initialize the dictionary
    sequencesByLength: Dict[int,List[str]] = dict()
    for readSize in readSizeRange: sequencesByLength[readSize] = list()

    # Populate the dictionary with sequences of appropriate lenth.
    with open(fastaFilePath, 'r') as fastaFile:
        for fastaEntry in FastaFileIterator(fastaFile):

            sequenceLength = len(fastaEntry.sequence)
            if sequenceLength in readSizeRange:
                sequencesByLength[sequenceLength].append(fastaEntry.sequence)

    return sequencesByLength


# Given one or more fasta files and the features to count, find which indices are enriched from the end of the sequences, by length.
def findEnrichedIndices(bedFilePaths: List[str], genomeFastaFilePath, fastaFilePaths: List[str],
                        countIndividualBases, countDipys, getSecondPlace):

    # Set some default parameters (I may change this to be more malleable at some point.)
    readSizeRange = range(11,41)
    fromStartValues = range(0) # These should be 1-based!  (0 can't be used because fromStart and fromEnd values are differentiated by +/-)
    fromEndValues = range(1,11)

    # First, convert any bed file paths to fasta format.
    print("Converting bed files to fasta format...")
    newFastaFilePaths = list()
    for bedFilePath in bedFilePaths:

        print(f"Converting {os.path.basename(bedFilePath)}...")
        intermediateDir = os.path.join(os.path.dirname(bedFilePath),"intermediate_files")
        checkDirs(intermediateDir)
        fastaBasename = os.path.basename(bedFilePath).rsplit('.',1)[0] + ".fa"
        fastaOutputFilePath = os.path.join(intermediateDir,fastaBasename)
        newFastaFilePaths.append(fastaOutputFilePath)

        bedToFasta(bedFilePath, genomeFastaFilePath, fastaOutputFilePath)
    
    # Add any new fasta file paths to the current list and then use them to search for enriched indices.
    fastaFilePaths += newFastaFilePaths
    for fastaFilePath in fastaFilePaths:

        print()
        print("Working with:",os.path.basename(fastaFilePath))

        # Generate a file path for the output file.
        if getIsolatedParentDir(fastaFilePath) == "intermediate_files":
            outputDir = os.path.dirname(os.path.dirname(fastaFilePath))
        else: outputDir = os.path.dirname(fastaFilePath)
        outputBasename = os.path.basename(fastaFilePath).rsplit('.',1)[0] + "_enriched_indices.tsv"
        enrichedIndicesOutputFilePath = os.path.join(outputDir,outputBasename)

        # Get a dictionary of sequences with lengths within readSizeRange
        print("Reading in sequences and binning by length...")
        sequencesByLength = getReadSequencesByLength(fastaFilePath, readSizeRange)

        # Next, prepare a dictionary to hold the final enriched indices info
        enrichedIndicesInfo: Dict[int,Dict[str,tuple]] = dict()
        if getSecondPlace: enrichedIndicesInfoSP: Dict[int,Dict[str,tuple]] = dict() # The "second place" dictionary.

        # For each sequence length, prepare a list of the enriched indices for each requested feature, as well as their frequency.
        print("Finding enriched indices for each sequence length bin...")
        for sequenceLength in readSizeRange:
            print("Working with sequences of length",sequenceLength)
            enrichedIndicesInfo[sequenceLength] = dict()
            if getSecondPlace: enrichedIndicesInfoSP[sequenceLength] = dict()

            if countIndividualBases:
                individualBaseFrequencyTable = BaseFrequencyTable(fromStartValues, fromEndValues, BaseFrequencyTable.TrackedFeature.singleBase)
                individualBaseFrequencyTable.generateBaseFrequencyTable(sequencesByLength[sequenceLength])
                for base in ('A','C','G','T'):
                    if getSecondPlace: enrichedIndicesInfoSP[sequenceLength][base] = individualBaseFrequencyTable.getMaxFrequencyAndPos(base, getSecondPlace)
                    enrichedIndicesInfo[sequenceLength][base] = individualBaseFrequencyTable.getMaxFrequencyAndPos(base)
            
            if countDipys:
                dipyFrequencyTable = BaseFrequencyTable(fromStartValues, fromEndValues, BaseFrequencyTable.TrackedFeature.dipys)
                dipyFrequencyTable.generateBaseFrequencyTable(sequencesByLength[sequenceLength])
                if getSecondPlace:
                    enrichedIndicesInfoSP[sequenceLength]["dipys"] = dipyFrequencyTable.getMaxFrequencyAndPos(BaseFrequencyTable.TrackedFeature.dipys, getSecondPlace)
                enrichedIndicesInfo[sequenceLength]["dipys"] = dipyFrequencyTable.getMaxFrequencyAndPos(BaseFrequencyTable.TrackedFeature.dipys)

        # Now, write the results to the output file!
        print("Writing Results...")
        with open(enrichedIndicesOutputFilePath, 'w') as enrichedIndicesOutputFile:

            # Write the header
            enrichedIndicesOutputFile.write("Sequence_Length" + '\t' + "Read_Count")
            for feature in enrichedIndicesInfo[readSizeRange[0]]:
                enrichedIndicesOutputFile.write('\t' + feature + "_Max_Frequency" + '\t' + feature + "_Max_Frequency_Position")
                if getSecondPlace: enrichedIndicesOutputFile.write('\t' + feature + "_Next_Max_Frequency" + '\t' + 
                                                                   feature + "_Next_Max_Frequency_Position" + '\t' + feature + "_Max_to_Next_Max_Diff")
            enrichedIndicesOutputFile.write('\n')

            # Write everything else!
            for sequenceLength in readSizeRange:
                enrichedIndicesOutputFile.write(str(sequenceLength) + '\t' + str(len(sequencesByLength[sequenceLength])))
                for feature in enrichedIndicesInfo[sequenceLength]:
                    maxFrequencyInfo = enrichedIndicesInfo[sequenceLength][feature]
                    enrichedIndicesOutputFile.write('\t' + str(maxFrequencyInfo[0]) + '\t' + str(maxFrequencyInfo[1]))

                    if getSecondPlace: 
                        nextMaxFrequencyInfo = enrichedIndicesInfoSP[sequenceLength][feature]
                        enrichedIndicesOutputFile.write('\t' + str(nextMaxFrequencyInfo[0]) + '\t' + str(nextMaxFrequencyInfo[1]) +
                                                        '\t' + str(maxFrequencyInfo[0] - nextMaxFrequencyInfo[0]))
                enrichedIndicesOutputFile.write('\n')


def main():

    # Create a simple dialog for selecting the fasta files and the nucleotide feature(s) to count.
    with TkinterDialog(workingDirectory=getDataDirectory()) as dialog:
        dialog.createMultipleFileSelector("Bed files of aligned data:", 0, "aligned_reads.bed", 
                                          ("Bed Files", ".bed"))
        dialog.createFileSelector("Genome fasta file (If bed files are provided)", 1,
                                  ("Fasta File", ".fa"))
        dialog.createMultipleFileSelector("Fasta files of aligned data:", 2, "aligned_reads.fa", 
                                          ("Fasta Files", ".fa"))
        dialog.createCheckbox("Count individual bases", 3, 0)
        dialog.createCheckbox("Count dipys", 3, 1)
        dialog.createCheckbox("Record second-most enriched positions", 4, 0)

    # Get the input for the findEnrichedIndices function
    bedFilePaths = dialog.selections.getFilePathGroups()[0]
    genomeFastaFilePath = dialog.selections.getIndividualFilePaths()[0]
    fastaFilePaths = dialog.selections.getFilePathGroups()[1]

    countIndividualBases = dialog.selections.getToggleStates()[0]
    countDipys = dialog.selections.getToggleStates()[1]
    getSecondPlace = dialog.selections.getToggleStates()[2]

    findEnrichedIndices(bedFilePaths, genomeFastaFilePath, fastaFilePaths,
                        countIndividualBases, countDipys, getSecondPlace)


if __name__ == "__main__": main()