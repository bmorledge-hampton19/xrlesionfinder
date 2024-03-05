import os
from enum import Enum
from typing import Dict, List
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
from benbiohelpers.FileSystemHandling.FastaFileIterator import FastaFileIterator
from benbiohelpers.FileSystemHandling.DirectoryHandling import checkDirs, getIsolatedParentDir, getTempDir
from benbiohelpers.FileSystemHandling.BedToFasta import bedToFasta
from xrlesionfinder.ProjectManagement.UsefulFileSystemFunctions import getDataDirectory


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


def findEnrichedIndices(bedFilePaths: List[str], genomeFastaFilePath, fastaFilePaths: List[str],
                        countIndividualBases, countDipys, getSecondPlace,
                        readSizeRange = range(16,36), fromStartValues = range(0,0), fromEndValues = range(1,16),
                        outputBulkFrequencies = False):
    """
    Given one or more fasta files and the features to count, find which indices are enriched for each sequence length.
    Right now, the search is restricted to specific read size ranges and positions relative to the sequence start and end:
    - readSizeRange specifies the reads that will be analyzed.
    - fromStartValues specifies the 1-based positions from the 5' end that will be analyzed.
    - fromEndValues specifies the 1-based positions from the 3- end that will be analyzed.
    - Default values reflect reasonable contraints for human XR-seq reads.

    Unfortunately, the code is pretty brittle at the moment. (e.g., depending on the above values, it may try to look up string indices that do not exist.)
    It is also very memory inefficient, due to reading all the valid fasta sequences into a dictionary at runtime... but it does the job!
    If I ever need to scale up this kind of analysis though, it may be better to redesign the core algorithm...
    """

    # Create a list of all searched positions.
    allSearchValues = list(fromStartValues) + ([-value for value in fromEndValues][::-1])

    # First, convert any bed file paths to fasta format.
    print("Converting bed files to fasta format...")
    newFastaFilePaths = list()
    for bedFilePath in bedFilePaths:

        print(f"Converting {os.path.basename(bedFilePath)}...")
        tmpDir = getTempDir(bedFilePath)
        checkDirs(tmpDir)
        fastaBasename = os.path.basename(bedFilePath).rsplit('.',1)[0] + ".fa"
        fastaOutputFilePath = os.path.join(tmpDir,fastaBasename)
        newFastaFilePaths.append(fastaOutputFilePath)

        bedToFasta(bedFilePath, genomeFastaFilePath, fastaOutputFilePath)
    
    # Add any new fasta file paths to the current list and then use them to search for enriched indices.
    fastaFilePaths += newFastaFilePaths
    for fastaFilePath in fastaFilePaths:

        print()
        print("Working with:",os.path.basename(fastaFilePath))

        # Generate output file paths.
        if getIsolatedParentDir(fastaFilePath) == ".tmp":
            outputDir = os.path.dirname(os.path.dirname(fastaFilePath))
        else: outputDir = os.path.dirname(fastaFilePath)
        outputBasename = os.path.basename(fastaFilePath).rsplit('.',1)[0] + "_enriched_indices.tsv"
        enrichedIndicesOutputFilePath = os.path.join(outputDir,outputBasename)

        if outputBulkFrequencies:
            if countIndividualBases:
                individualFrequenciesOutputFilePath = os.path.join(outputDir, os.path.basename(fastaFilePath).rsplit('.',1)[0] + "_individual_nuc_frequencies.tsv")
                individualFrequenciesOutputFile = open(individualFrequenciesOutputFilePath, 'w')
                individualFrequenciesOutputFile.write('\t'.join(("Sequence_Length", "Position", "A_Frequency", "C_Frequency", "G_Frequency", "T_Frequency")) + '\n')
            if countDipys:
                dipyFrequenciesOutputFilePath = os.path.join(outputDir, os.path.basename(fastaFilePath).rsplit('.',1)[0] + "_dipy_frequencies.tsv")
                dipyFrequenciesOutputFile = open(dipyFrequenciesOutputFilePath, 'w')
                dipyFrequenciesOutputFile.write('\t'.join(("Sequence_Length", "Position", "CC_Frequency", "CT_Frequency", "TC_Frequency", "TT_Frequency")) + '\n')

        # Get a dictionary of sequences with lengths within readSizeRange. (This is horribly memory inefficient...)
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
                
                # If requested, write the individual frequencies.
                if outputBulkFrequencies:
                    for searchValue in allSearchValues:
                        individualFrequenciesOutputFile.write('\t'.join((str(sequenceLength), str(searchValue),
                                                                    str(individualBaseFrequencyTable.baseFrequencies['A'][searchValue]),
                                                                    str(individualBaseFrequencyTable.baseFrequencies['C'][searchValue]),
                                                                    str(individualBaseFrequencyTable.baseFrequencies['G'][searchValue]),
                                                                    str(individualBaseFrequencyTable.baseFrequencies['T'][searchValue]),
                        )) + '\n')
                
            
            if countDipys:
                dipyFrequencyTable = BaseFrequencyTable(fromStartValues, fromEndValues, BaseFrequencyTable.TrackedFeature.dipys)
                dipyFrequencyTable.generateBaseFrequencyTable(sequencesByLength[sequenceLength])
                if getSecondPlace:
                    enrichedIndicesInfoSP[sequenceLength]["dipys"] = dipyFrequencyTable.getMaxFrequencyAndPos(BaseFrequencyTable.TrackedFeature.dipys, getSecondPlace)
                enrichedIndicesInfo[sequenceLength]["dipys"] = dipyFrequencyTable.getMaxFrequencyAndPos(BaseFrequencyTable.TrackedFeature.dipys)

                # If requested, write the individual counts.
                if outputBulkFrequencies:
                    for searchValue in allSearchValues:
                        individualFrequenciesOutputFile.write('\t'.join((str(sequenceLength), str(searchValue),
                                                                    str(individualBaseFrequencyTable.baseFrequencies['CC'][searchValue]),
                                                                    str(individualBaseFrequencyTable.baseFrequencies['CT'][searchValue]),
                                                                    str(individualBaseFrequencyTable.baseFrequencies['TC'][searchValue]),
                                                                    str(individualBaseFrequencyTable.baseFrequencies['TT'][searchValue]),
                        )) + '\n')

        # Close bulk frequency ouptut file paths as necessary.
        if outputBulkFrequencies and countIndividualBases: individualFrequenciesOutputFile.close()
        if outputBulkFrequencies and countDipys: dipyFrequenciesOutputFile.close()

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

    with TkinterDialog(workingDirectory = getDataDirectory(), title = "Find Enriched Indices") as dialog:
        dialog.createMultipleFileSelector("Bed files of aligned data:", 0, "aligned_reads.bed", 
                                          ("Bed Files", ".bed"))
        dialog.createGenomeSelector(1, 0)
        dialog.createMultipleFileSelector("Fasta files of aligned data:", 2, "aligned_reads.fa", 
                                          ("Fasta Files", ".fa"))
        dialog.createCheckbox("Count individual bases", 3, 0)
        dialog.createCheckbox("Count dipys", 3, 1)
        dialog.createCheckbox("Record second-most enriched positions", 4, 0)
        dialog.createCheckbox("Output bulk frequencies", 4, 1)

    # Get the input for the findEnrichedIndices function
    bedFilePaths = dialog.selections.getFilePathGroups()[0]
    genomeFastaFilePath = dialog.selections.getGenomes(returnType="fasta")[0]
    fastaFilePaths = dialog.selections.getFilePathGroups()[1]

    countIndividualBases = dialog.selections.getToggleStates()[0]
    countDipys = dialog.selections.getToggleStates()[1]
    getSecondPlace = dialog.selections.getToggleStates()[2]
    outputBulkFrequencies = dialog.selections.getToggleStates()[3]

    findEnrichedIndices(bedFilePaths, genomeFastaFilePath, fastaFilePaths,
                        countIndividualBases, countDipys, getSecondPlace,
                        outputBulkFrequencies = outputBulkFrequencies)


if __name__ == "__main__": main()