# This script will be called from the command line to execute other scripts.
from argparse import ArgumentParser
from benbiohelpers.CustomErrors import *
from xrlesionfinder.AlignmentAndFormatting import AlignXRSeqReads
from xrlesionfinder.ProjectManagement.GenomeManager import GenomeManagerError
import argparse, importlib.util, sys, traceback
if importlib.util.find_spec("shtab") is not None: 
        import shtab
        fileCompletion = shtab.FILE
else:
        fileCompletion = None


def formatAlignReadsParser(alignReadsParser: ArgumentParser):

    alignReadsParser.set_defaults(func = AlignXRSeqReads.parseArgs)
    alignReadsParser.add_argument("XRSeqReadsFilePaths", nargs = '*',
                                  help = "One or more paths to XR-seq reads files (in fastq format). Can be gzipped. "
                                         "If given a directory, it will be recursively searched for files ending in \".fastq\" or "
                                         "\".fastq.gz\".").complete = fileCompletion
    # TODO: Finish this.


def getMainParser():

    # Initialize the argument parser.
    parser = argparse.ArgumentParser(description = "Run xrlesionfinder on XR-seq reads to identify lesion positions.",
                                     prog = "xrlesionfinder")
    subparsers = parser.add_subparsers(title='xrlesionfinder command', required = True, dest = "xrlesionfinder command")

    ### Create the subparsers for each relevant script.

    # For obtaining reads...
    # TODO: Make parser for this ^

    # For aligning reads...
    alignReadsParser = subparsers.add_parser("alignreads", description = "Align XR-seq reads in preparation for identifying lesions.")
    formatAlignReadsParser(alignReadsParser)
    

    return parser


def main():

    # Run the relevant function for the subparser given.
    args = getMainParser().parse_args()
    try:
        args.func(args)
    except MetadataPathError as error:
        sys.exit("Error finding metadata expected at:\n" + error.path + "\n"
                 "Make sure that the related directory was created through mutperiod and that "
                 "you have not manually altered the file structure within the \"mutperiod_data\" directory.")
    except UserInputError as error:
        sys.exit("Error: " + str(error))
    except GenomeManagerError as error:
        sys.exit(f"Error: {error}\n Use the command \"xrlesionfinder addgenome\" to add/update genome locations.")
    except Exception:
        traceback.print_exc()
        print("\n\n\n")
        sys.exit("Unexpected error encountered.  For more assistance, please send the above traceback along with "
                 "an explanation of what caused the error to b.morledge-hampton@wsu.edu")



if __name__ == "__main__": main()