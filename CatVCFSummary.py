#!/usr/bin/env python2.7
import os
import argparse
import logging
import sys

from Utils import configure_logging
from SummarizeVCF import VCFSummaryParser, merge_two_vcf_summaries

def configure_argparser(argparser_obj):

    def file_type(arg_string):
        """
        This function check both the existance of input file and the file size
        :param arg_string: file name as string
        :return: file name as string
        """
        if not os.path.exists(arg_string):
            err_msg = "%s does not exist! " \
                      "Please provide a valid file!" % arg_string
            raise argparse.ArgumentTypeError(err_msg)

        return arg_string

    # Path to VCF input file
    argparser_obj.add_argument("-i",
                               action="store",
                               type=file_type,
                               nargs='+',
                               dest="input_files",
                               required=True,
                               help="Space-delimited list of RecodedVCF files to combine")

    # Verbosity level
    argparser_obj.add_argument("-v",
                               action='count',
                               dest='verbosity_level',
                               required=False,
                               default=0,
                               help="Increase verbosity of the program."
                                    "Multiple -v's increase the verbosity level:\n"
                                    "0 = Errors\n"
                                    "1 = Errors + Warnings\n"
                                    "2 = Errors + Warnings + Info\n"
                                    "3 = Errors + Warnings + Info + Debug")

def main():

    # Configure argparser
    argparser = argparse.ArgumentParser(prog="CatVCFSummary")
    configure_argparser(argparser)

    # Parse the arguments
    args = argparser.parse_args()

    # Configure logging
    configure_logging(args.verbosity_level)

    try:

        # Create VCFSummaryParser objects from files
        vcf_summary_parsers = [VCFSummaryParser(input_file) for input_file in args.input_files]

        # Parse and merge VCFSummaries
        merged_summary = None
        for vcf_summary_parser in vcf_summary_parsers:

            # Parse VCFSummary
            vcf_summary = vcf_summary_parser.parse()

            if merged_summary is None:
                # Initialize base VCFSummary
                merged_summary = vcf_summary
            else:
                # Merge with previous VCFSummaries
                merged_summary = merge_two_vcf_summaries(merged_summary, vcf_summary)

        print merged_summary

    except KeyboardInterrupt:
        logging.error("(Main) Keyboard interrupt!")
        raise

    except BaseException, e:
        # Report any errors that arise
        logging.error("(Main) CatVCFSummary failed!")
        if e.message != "":
            logging.error("Received the following error message:\n%s" % e.message)
        raise

if __name__ == "__main__":
    sys.exit(main())
