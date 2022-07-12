# This script defines the metadata structure for the xrlesionfinder package.
import os
from benbiohelpers.DataPipelineManagement.Metadata import Metadata, MetadataFeatureID, MetadataFeatureValue
from enum import auto

class DataType(MetadataFeatureValue):
    """
    Defines what kind of data the file holds.
    Values are used to enforce file naming conventions.
    Each value is a tuple of the data type string and the file extension.
    """
    READS_FASTQ_FILE = auto(), ("",".fastq")
    ALIGNMENT_SAM_FILE = auto(), ("",".sam")
    ALIGNMENT_BAM_FILE = auto(), ("",".bam")
    ALIGNMENT_BED_FILE = auto(), ("",".bed")


class XRLFMetadataFeatureID(MetadataFeatureID):
    ORGANISM = auto(), str
    LESION = auto(), str
    TIMEPOINT = auto(), str
    REPETITION = auto(), str
    DATA_TYPE = auto(), DataType
XRLFMFID = XRLFMetadataFeatureID

class XRLFMetadata(Metadata):
    FeatureIDEnum = XRLFMetadataFeatureID

    defaultValues = {
        XRLFMFID.REPETITION:""
    }

    def getFilePath(self, useParentDirectory = True):

        filePathPieces = [self[XRLFMFID.ORGANISM], self[XRLFMFID.LESION],
                          self[XRLFMFID.TIMEPOINT], self[XRLFMFID.REPETITION]]

        # Remove any empty string file path pieces and make sure we don't have any nonetypes.
        assert None not in filePathPieces, "NoneType in filePathPieces. Is something missing a default value?"
        filePathPieces = [str(filePathPiece) for filePathPiece in filePathPieces if str(filePathPiece)]

        if useParentDirectory: directory = os.path.dirname(self.directory)
        else: directory = self.directory

        return os.path.join(directory,'_'.join(filePathPieces) + ''.join(self[XRLFMFID.DATA_TYPE].value))