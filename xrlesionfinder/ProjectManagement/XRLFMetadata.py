# This script defines the metadata structure for the xrlesionfinder package.
import os
from typing import List
from benbiohelpers.DataPipelineManagement.Metadata import Metadata, MetadataFeatureID, MetadataFeatureValue
from enum import auto

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
    ALT_ID = auto(), str
    ORGANISM = auto(), str
    LESION = auto(), str
    TIMEPOINT = auto(), str
    REPETITION = auto(), str
    DATA_TYPE = auto(), DataType
    MISMATCHES = auto(), list
    EXPANSION_NUM = auto(), int
    FILTERING = auto(), str
XRLFMFID = XRLFMetadataFeatureID

class XRLFMetadata(Metadata):
    FeatureIDEnum = XRLFMetadataFeatureID

    defaultValues = {
        XRLFMFID.MISMATCHES: list(),
        XRLFMFID.EXPANSION_NUM: 0,
        XRLFMFID.FILTERING: '',
    }

    def getFilePath(self, useParentDirectory = True):

        # First, check for an alternate ID. If it's found, it means that
        # the individual components that usually identify the data set couldn't be found
        # and the alternate ID should be used in their place.
        if self[XRLFMFID.ALT_ID] is not None:
            filePathPieces = [self[XRLFMFID.ALT_ID]]
        else:
            filePathPieces = [self[XRLFMFID.ORGANISM], self[XRLFMFID.LESION],
                              self[XRLFMFID.TIMEPOINT], self[XRLFMFID.REPETITION]]

        # Incoroporate information on mismatches up to 2; afterwards, the string just becomes "multi_mismatch".
        mismatches: List[str] = self[XRLFMFID.MISMATCHES]
        if len(mismatches < 3): filePathPieces += [mismatch.replace('>', "_to_") for mismatch in mismatches]
        else: filePathPieces.append("multi_mismatch")

        if self[XRLFMFID.EXPANSION_NUM] > 0: filePathPieces.append(f"{self[XRLFMFID.EXPANSION_NUM]}bp_expanded")

        filePathPieces.append(self[XRLFMFID.FILTERING])

        # Remove any empty string file path pieces and make sure we don't have any nonetypes.
        assert None not in filePathPieces, "NoneType in filePathPieces. Is something missing a default value?"
        filePathPieces = [str(filePathPiece) for filePathPiece in filePathPieces if str(filePathPiece)]

        if useParentDirectory: directory = os.path.dirname(self.directory)
        else: directory = self.directory

        return os.path.join(directory,'_'.join(filePathPieces) + '_' + ''.join(self[XRLFMFID.DATA_TYPE].value))


