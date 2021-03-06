#!/usr/bin/env python2.7
import os
import argparse
import logging
import sys

from VCF import VCFHelper
from Utils import configure_logging
from SummarizeVCF import VCFSummarizer
from SummarizeVCF.VariantAnalyzer import VariantAnalyzerFactory

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

    # Type of VCFSummary to get
    argparser_obj.add_argument('summary_type',
                               choices=VariantAnalyzerFactory.SUMMARY_TYPES.keys(),
                               help="VCF Summary type.")

    # Path to VCF input file
    argparser_obj.add_argument("--vcf",
                               action="store",
                               type=file_type,
                               dest="vcf_file",
                               required=True,
                               help="Path to vcf file to recode.")

    # Upper boundary of indel length summary
    argparser_obj.add_argument("--max-records",
                               action="store",
                               type=int,
                               dest="max_records",
                               required=False,
                               default=-1,
                               help="Maximum number of records to process. Default: ALL.")

    # Upper boundary of indel length summary
    argparser_obj.add_argument("--max-indel-len",
                               action="store",
                               type=int,
                               dest="max_indel_len",
                               required=False,
                               default=100,
                               help="Upper bound of indel length summary.")

    # Upper boundary of indel length summary
    argparser_obj.add_argument("--max-depth",
                               action="store",
                               type=int,
                               dest="max_depth",
                               required=False,
                               default=100,
                               help="Upper bound of variant depth summary.")

    # Upper boundary of indel length summary
    argparser_obj.add_argument("--max-qual",
                               action="store",
                               type=int,
                               dest="max_qual",
                               required=False,
                               default=100,
                               help="Upper bound of variant quality summary.")

    # Number of bins for allele frequency spectrum
    argparser_obj.add_argument("--afs-bins",
                               action="store",
                               type=int,
                               dest="num_afs_bins",
                               required=False,
                               default=20,
                               help="Number of bins to use for Allele Frequency Spectrum.")

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
    argparser = argparse.ArgumentParser(prog="SummarizeVCF",
                                        usage="SummarizeVCF.py <summary_type> [options]")
    configure_argparser(argparser)

    # Parse the arguments
    args = argparser.parse_args()

    # Configure logging
    configure_logging(args.verbosity_level)

    # Get names of input/output files
    vcf_file                = args.vcf_file
    summary_type            = args.summary_type
    max_records             = args.max_records

    # Args for the VCF summary
    summary_args = {"max_indel_len" : args.max_indel_len}
    summary_args["max_depth"]     = args.max_depth
    summary_args["max_qual"]      = args.max_qual
    summary_args["num_afs_bins"]  = args.num_afs_bins

    try:

        logging.debug("(Main) Starting to summarize VCF file: %s" % vcf_file)

        # Check to make sure VCF file is valid
        if not VCFHelper.is_valid_vcf(vcf_file):
            raise IOError("Invalid VCF file!")

        # Summarize VCF and print to stdout
        summarizer = VCFSummarizer(vcf_file, summary_type, max_records, **summary_args)
        summarizer.summarize()
        print summarizer.get_summary()

        # Summarize VCF file and print output to outfile
        logging.debug("(Main) Successfully summarized VCF file!")

    except KeyboardInterrupt:
        logging.error("(Main) Keyboard interrupt!")
        raise

    except BaseException, e:
        # Report any errors that arise
        logging.error("(Main) VcfQC failed!")
        if e.message != "":
            logging.error("Received the following error message:\n%s" % e.message)
        raise

if __name__ == "__main__":
    sys.exit(main())
