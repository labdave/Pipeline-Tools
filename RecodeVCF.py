#!/usr/bin/env python2.7
import os
import argparse
import logging
import sys

from FileValidators import VCFValidator
from Utils import configure_logging

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
    argparser_obj.add_argument("--vcf",
                               action="store",
                               type=file_type,
                               dest="vcf_file",
                               required=True,
                               help="Path to vcf file to recode.")

    # Path to recoded output file
    argparser_obj.add_argument("--min-call-depth",
                               action="store",
                               type=int,
                               dest="min_call_depth",
                               required=True,
                               default=10,
                               help="Minimum read depth required to call variant genotype.")

    # Path to recoded output file
    argparser_obj.add_argument("--missing-data-char",
                               action="store",
                               type=str,
                               dest="missing_data_char",
                               required=True,
                               default='.',
                               help="Minimum read depth required to call variant genotype.")

    # Path to recoded output file
    argparser_obj.add_argument("--output",
                               action="store",
                               type=str,
                               dest="output_file",
                               required=True,
                               help="Path to output file.")

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
    argparser = argparse.ArgumentParser(prog="RecodeVCF")
    configure_argparser(argparser)

    # Parse the arguments
    args = argparser.parse_args()

    # Configure logging
    configure_logging(args.verbosity_level)

    # Get names of input/output files
    vcf_file            = args.vcf_file
    out_file            = args.out_file
    min_call_depth      = args.min_call_depth
    missing_data_char   = args.missing_data_char

    try:

        logging.debug("(Main) Starting to recode VCF file: %s" % vcf_file)

        # Check to make sure VCF file is valid
        if not VCFValidator.is_valid(vcf_file):
            raise IOError("Invalid VCF file!")

        # Stuff goes here
        logging.debug("(Main) Successfully recoded VCF file to output: %s" % out_file)

    except KeyboardInterrupt:
        logging.error("(Main) Keyboard interrupt!")
        raise

    except BaseException, e:
        # Report any errors that arise
        logging.error("(Main) RecodeVCF failed!")
        if e.message != "":
            logging.error("Received the following error message:\n%s" % e.message)
        raise

if __name__ == "__main__":
    sys.exit(main())
