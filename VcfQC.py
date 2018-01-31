#!/usr/bin/env python2.7
import os
import argparse
import logging
import sys
import vcf

from VCF import VCFHelper
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

    # Path to VCF input file
    argparser_obj.add_argument("--output",
                               action="store",
                               type=str,
                               dest="out_file",
                               required=True,
                               help="Path to recoded output file.")

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
                               default=500,
                               help="Upper bound of variant depth summary.")

    # Upper boundary of indel length summary
    argparser_obj.add_argument("--max-qual",
                               action="store",
                               type=int,
                               dest="max_qual",
                               required=False,
                               default=250,
                               help="Upper bound of variant quality summary.")

    # Number of bins for allele frequency spectrum
    argparser_obj.add_argument("--afs-bins",
                               action="store",
                               type=int,
                               dest="num_afs_bins",
                               required=False,
                               default=20,
                               help="Number of bins to use for Allele Frequency Spectrum.")

    # Character used to denote missing variant information
    argparser_obj.add_argument("--missing-data-char",
                               action="store",
                               type=str,
                               dest="missing_data_char",
                               required=False,
                               default='.',
                               help="Character used as placeholder for missing VCF info.")

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

def summarize_vcf(vcf_parser, vcf_summary):

    # Line counter
    processed = 0

    for record in vcf_parser:

        if processed % 100000 == 0:
            logging.info("Processed %d records!" % processed)

        # Check to make sure there is only one alternate allele for the record
        if len(record.ALT) > 1:
            # Raise error because it's just easier if we make people normalize their data beforehand
            logging.error("Multiallelic error! All VCF records must have only one alternate allele. "
                          "Received the following record:\n%s" % record)
            # Failing to raise this error could cause some weirdo unintended bugs
            raise IOError("Multiple alternate alleles detected at one or more positions in vcf file!")

        # Process variant record and add pertinant data to summary
        vcf_summary.process_record(record)

def main():

    # Configure argparser
    argparser = argparse.ArgumentParser(prog="VcfQC")
    configure_argparser(argparser)

    # Parse the arguments
    args = argparser.parse_args()

    # Configure logging
    configure_logging(args.verbosity_level)

    # Get names of input/output files
    vcf_file                = args.vcf_file
    out_file                = args.out_file
    max_indel_len           = args.max_indel_len
    max_depth               = args.max_depth
    max_qual                = args.max_qual
    num_afs_bins            = args.num_afs_bins
    missing_data_char       = args.missing_data_char

    try:

        logging.debug("(Main) Starting to summarize VCF file: %s" % vcf_file)

        # Check to make sure VCF file is valid
        if not VCFHelper.is_valid_vcf(vcf_file):
            raise IOError("Invalid VCF file!")

        # Get VCFParser
        logging.debug("Initializing VCF parser for file: %s" % vcf_file)
        vcf_parser = vcf.Reader(open(vcf_file, "r"))

        # Initialize VCF summary
        #sample_names    = vcf_parser.samples

        #info_fields = OrderedDict.fromkeys(vcf_parser.infos.keys()).keys()
        info_fields = VCFHelper.get_info_field_names(vcf_parser)
        logging.debug("Available annotation fields:\n%s" % info_fields)

        #vcf_summary     = VCFSummary(sample_names, max_depth, max_qual, max_indel_len, num_afs_bins)
        #logging.debug("Samples:\n%s" % ", ".join(sample_names))

        # Parse VCF and create variant summary
        #summarize_vcf(vcf_parser, vcf_summary)

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
