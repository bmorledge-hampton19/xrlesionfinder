# This script defines the metadata structure for the xrlesionfinder package.
import os, warnings
from typing import List
from benbiohelpers.DataPipelineManagement.Metadata import Metadata, MetadataFeatureID, MetadataFeatureValue
from benbiohelpers.CustomErrors import MetadataAutoGenerationError
from enum import auto

def willItFloat(floatContestant):
    try: float(floatContestant)
    except ValueError: return False
    else: return True

expectedCellTypes = ["NHF1", "GM12878", "yeast", "Dm", "Ecoli"]
expectedLesions = ["CPD", "6-4", "6-4PP", "cisplatin"]

class DataType(MetadataFeatureValue):
    """
    Defines what kind of data the file holds.
    Values are used to enforce file naming conventions.
    Each value is a tuple of the data type string and the file extension.
    """
    READS_FASTQ_FILE = auto(), ('',".fastq")
    TRIMMED_READS_FASTQ_FILE = auto(), ("trimmed",".fastq")
    ALIGNMENT_SAM_FILE = auto(), ('',".sam")
    ALIGNMENT_BAM_FILE = auto(), ('',".bam")
    ALIGNMENT_BED_FILE = auto(), ('',".bed")

class XRLFMetadataFeatureID(MetadataFeatureID):
    ALT_ID = auto(), str # e.g. "weird_organism_mystery_lesion_23_days_rrepp_1.5"
    CELL_TYPE = auto(), str # e.g. "NHF1" or "Drosophila_melanogaster"
    LESION = auto(), str # e.g. "CPD" or "6-4"
    TIMEPOINT = auto(), str # e.g. "10m" or "2h"
    REPETITION = auto(), str # e.g. "rep1" or "all_reps"
    MISMATCHES = auto(), list # e.g. ["C>T", "C>G", "C>A"]
    EXPANSION_NUM = auto(), int # e.g. 10
    STRAND_POLARITY = auto(), str # three_prime or five_prime
    FILTERING = auto(), str # e.g. "TGG_filtered"
    DATA_TYPE = auto(), DataType # See enum above
XRLFMFID = XRLFMetadataFeatureID

class XRLFMetadata(Metadata):
    FeatureIDEnum = XRLFMetadataFeatureID

    defaultValues = {
        XRLFMFID.MISMATCHES: list(),
        XRLFMFID.EXPANSION_NUM: 0,
        XRLFMFID.FILTERING: '',
        XRLFMFID.STRAND_POLARITY: '',
    }

    def getFilePath(self, useParentDirectory = True):

        # First, check for an alternate ID. If it's found, it means that
        # the individual components that usually identify the data set couldn't be found
        # and the alternate ID should be used in their place.
        if self[XRLFMFID.ALT_ID] is not None:
            filePathPieces = [self[XRLFMFID.ALT_ID]]
        else:
            filePathPieces = [self[XRLFMFID.CELL_TYPE], self[XRLFMFID.LESION],
                              self[XRLFMFID.TIMEPOINT], self[XRLFMFID.REPETITION]]

        # Incoroporate information on mismatches up to 2; afterwards, the string just becomes "multi_mismatch".
        mismatches: List[str] = self[XRLFMFID.MISMATCHES]
        if len(mismatches < 3): filePathPieces += [mismatch.replace('>', "_to_") for mismatch in mismatches]
        else: filePathPieces.append("multi_mismatch")

        if self[XRLFMFID.EXPANSION_NUM] > 0: filePathPieces.append(f"{self[XRLFMFID.EXPANSION_NUM]}bp_expanded")

        filePathPieces.append(self[XRLFMFID.FILTERING])
        filePathPieces.append(self[XRLFMFID.STRAND_POLARITY])

        # Remove any empty string file path pieces and make sure we don't have any nonetypes.
        assert None not in filePathPieces, "NoneType in filePathPieces. Is something missing a default value?"
        filePathPieces = [str(filePathPiece) for filePathPiece in filePathPieces if str(filePathPiece)]

        if useParentDirectory: directory = os.path.dirname(self.directory)
        else: directory = self.directory

        dataTypeString = ''.join(self[XRLFMFID.DATA_TYPE].value)
        if not dataTypeString.startswith('.'): dataTypeString = '_' + dataTypeString

        return os.path.join(directory,'_'.join(filePathPieces) + dataTypeString)


    def getFeaturesFromString(self, featuresString: str):
        
        # Eventually this function should attempt to derive an organism, lesion, timepoint, and repetition ID from the string.
        # These features should be at the beginning but are generally more free-form, so it is easiest if all other parts of the
        # string are parsed and removed first. To this end, parsing will start at the end and work backwards.

        # First, the data type needs to be determined. data types can contain one another as substrings, so
        # all data types will need to be tested with "endswith" and the longest will be used.
        longestDataType = None
        dataTypeStringLength = 0
        for dataType in DataType:
            dataTypeString = ''.join(dataType.value)
            if featuresString.endswith(dataTypeString):
                if len(dataTypeString) > dataTypeStringLength:
                    dataTypeStringLength = len(dataTypeString)
                    longestDataType = dataType
        if longestDataType is None:
            raise MetadataAutoGenerationError(f"Could not find a valid data type on the end of {featuresString}.")
        else:
            self[XRLFMFID.DATA_TYPE] = longestDataType
            featuresString = featuresString.rsplit(''.join(longestDataType.value), 1)[0]
            if featuresString.endswith('_'): featuresString = featuresString[:-1]

        # Next, look for optional information: Expansion number, filtering, and mismatch info.
        if featuresString.endswith("three_prime"):
            self[XRLFMFID.STRAND_POLARITY] = "three_prime"
            featuresString = featuresString.rsplit("_three_prime", 1)[0]
        elif featuresString.endswith("five_prime"):
            self[XRLFMFID.STRAND_POLARITY] = "five_prime"
            featuresString = featuresString.rsplit("_five_prime", 1)[0]

        if featuresString.endswith("filtered"):
            filtering = featuresString.split('_')[-2] + "_filtered"
            self[XRLFMFID.FILTERING] = filtering
            featuresString = featuresString.rsplit('_' + filtering, 1)[0]

        if featuresString.endswith("_bp_expanded"):
            expansionNum = int(featuresString.split('_')[-3])
            self[XRLFMFID.EXPANSION_NUM] = expansionNum
            featuresString = featuresString.rsplit(f"_{expansionNum}_bp_expanded", 1)[0]

        if featuresString.endswith("multi_mismatch"):
            self[XRLFMFID.MISMATCHES] = ["multi_mismatch"]
            featuresString = featuresString.rsplit("_multi_mismatch", 1)[0]
        else:
            splitFeaturesString = featuresString.split('_')
            mismatches = list()
            while (len(splitFeaturesString) >= 4 and splitFeaturesString[-2] == "to" and 
                   {'A','C','T','G'}.issuperset(set((splitFeaturesString[-1] + splitFeaturesString[-3]).upper()))):
                mismatches.append(f"{splitFeaturesString[-3]}>{splitFeaturesString[-1]}")
                featuresString = featuresString.rsplit('_',3)[0]
                splitFeaturesString = featuresString.split('_')
            self[XRLFMFID.MISMATCHES] = mismatches

        # Ok, NOW we can search for the organism, lesion, timepoint, and repetition. 
        # If this search fails at any point, everything we have left just gets lumped into the ALT_ID metadata
        # feature and throw a warning. (We tried!)
        potentialAltID = featuresString
        try:

            # First, search for the timepoint since it should actually have a somewhat regular pattern. (ends with 
            # 'm' or 'h' and beginning can be cast to a float.)
            potentialTimepoint = None
            for i,part in enumerate(featuresString.split('_')):
                if (part.lower().endswith('h') or part.lower().endswith('m') and willItFloat(part[:-1])):
                    if potentialTimepoint is not None: raise MetadataAutoGenerationError("Multiple potential timepoints found.")
                    elif i == 0:
                        raise MetadataAutoGenerationError("Timepoint found at beginning of ID string. "
                                                          "Where are the cell type and lesion IDs?")
                    elif i == featuresString.count('_'):
                        raise MetadataAutoGenerationError("Timepoint found at end of ID string. "
                                                          "Where is the repetition ID?")
                    else: potentialTimepoint = part
            if potentialTimepoint is None: raise MetadataAutoGenerationError("No timepoint information found.")
            self[XRLFMFID.TIMEPOINT] = potentialTimepoint

            # If we got this far, we have our timepoint! Split on it to get the repetition on one side and the
            # cell type/lesion on the other.
            featuresString, repetition = featuresString.split('_'+potentialTimepoint+'_')
            self[XRLFMFID.REPETITION] = repetition

            # Finally, look for expected lesions and cell types to decouple them.
            decoupled = False
            for lesion in expectedLesions:
                if featuresString.endswith(lesion):
                    self[XRLFMFID.LESION] = lesion
                    self[XRLFMFID.CELL_TYPE] = featuresString.rsplit('_' + lesion, 1)[0]
                    decoupled = True
                    break
            if not decoupled:
                for cellType in expectedCellTypes:
                    if featuresString.startswith(cellType):
                        self[XRLFMFID.CELL_TYPE] = cellType
                        self[XRLFMFID.LESION] = featuresString.split(cellType + '_', 1)[1]
                        decoupled = True
                        break
            
            if not decoupled: raise MetadataAutoGenerationError("Unable to find expected lesion or cell type to decouple them.")

        except MetadataAutoGenerationError as error:
            print(f"Unable to fully disambiguate metadata from given string: {error}\n Assigning remaining string to alt ID.")
            self[XRLFMFID.ALT_ID] = potentialAltID
