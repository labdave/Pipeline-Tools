#!/usr/bin/env python2.7
import os
import argparse
import logging
import sys

from VCF import VCFHelper, VCFAnnotationType, VCFRecoder, SnpEffVCFRecoder
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

    argparser_obj.add_argument("--info-column-file",
                               action="store",
                               type=file_type,
                               dest="info_column_file",
                               required=False,
                               default=None,
                               help="Path to file containing INFO columns to include in output (One per line.)")

    # Path to recoded output file
    argparser_obj.add_argument("--min-call-depth",
                               action="store",
                               type=int,
                               dest="min_call_depth",
                               required=False,
                               default=10,
                               help="Minimum read depth required to call variant genotype.")

    # Character used to denote missing variant information
    argparser_obj.add_argument("--missing-data-char",
                               action="store",
                               type=str,
                               dest="missing_data_char",
                               required=False,
                               default='.',
                               help="Character used as placeholder for missing VCF info.")

    # Character used to denote missing genotypes
    argparser_obj.add_argument("--missing-gt-char",
                               action="store",
                               type=str,
                               dest="missing_gt_char",
                               required=False,
                               default='NA',
                               help="Character used as placeholder for missing genotypes.")

    # Flag to allow VCF records with >1 alternate allele. Default is to not use this and throw errors.
    # You really shouldn't ever use this flag because parsing multiple alternate alleles with multiple annotations can fuck shit up.
    argparser_obj.add_argument("--multiallelic",
                               action="store_true",
                               dest="multiallelic",
                               required=False,
                               help="Flag allowing variant records to contain more than one alternate allele. This flag shouldn't really be used.")

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
    vcf_file                = args.vcf_file
    out_file                = args.out_file
    info_column_file        = args.info_column_file
    min_call_depth          = args.min_call_depth
    missing_data_char       = args.missing_data_char
    missing_gt_char         = args.missing_gt_char
    multiallelic            = args.multiallelic

    if multiallelic:
        # Throw warning about using the multiallelic mode
        logging.warning("**********DISCLAIMER**********\n"
                        "RecodeVCF is set to ALLOW sites with multiple alternate alleles. This can lead to strange behavior when parsing annotations.\n"
                        "We highly recommend you do NOT use this setting and instead use the bcftools norm function to split multiallelic sites into"
                        "multiple variant records!!!!\n")
    try:

        logging.debug("(Main) Starting to recode VCF file: %s" % vcf_file)

        # Check to make sure VCF file is valid
        if not VCFHelper.is_valid_vcf(vcf_file):
            raise IOError("Invalid VCF file!")

        # Get correct VCFparser based on annotation type
        vcf_recoder = None
        vcf_anno_type = VCFHelper.get_annotation_type(vcf_file)
        logging.debug("Detected VCF annotation type: '%s'" % vcf_anno_type)

        if vcf_anno_type == VCFAnnotationType.ANNOVAR or vcf_anno_type == VCFAnnotationType.UNKNOWN:
            # Get annovar/base VCFRecoder
            vcf_recoder = VCFRecoder(vcf_file,
                                     out_file,
                                     info_column_file = info_column_file,
                                     min_call_depth = min_call_depth,
                                     missing_data_char = missing_data_char,
                                     missing_gt_char = missing_gt_char,
                                     multiallelic=multiallelic)

        elif vcf_anno_type == VCFAnnotationType.SNPEFF:
            # Get snpeff VCFRecoder
            vcf_recoder = SnpEffVCFRecoder(vcf_file,
                                           out_file,
                                           info_column_file=info_column_file,
                                           min_call_depth = min_call_depth,
                                           missing_data_char = missing_data_char,
                                           missing_gt_char = missing_gt_char,
                                           multiallelic=multiallelic)

        # Recode the VCF file and print output to stdout
        vcf_recoder.recode_vcf()

        # Stuff goes here
        logging.debug("(Main) Successfully recoded VCF file!")
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
