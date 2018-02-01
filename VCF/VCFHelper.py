import logging
import vcf
from collections import OrderedDict

from VCFAnnotationType import VCFAnnotationType

class VCFHelper:
    # Collection of static functions relating to VCF files

    # Canoncial columns in the order they're supposed to appear in a VCF
    HEADER_FIELDS = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]

    @staticmethod
    def is_valid_vcf(path):

        num_samples = 0

        # Check to make sure VCF file is valid format
        with open(path, "r") as vcf_fh:

            # Boolean state variable determining which checks have been performed
            state = "FormatDeclaration"

            for line in vcf_fh:

                if state == "FormatDeclaration":

                    # Check to make sure the first line says ##fileformat=VCFv4.X
                    if not line.startswith("##fileformat=VCF"):
                        logging.error("Invalid VCF file format! First line of file doesn't look like '##fileformat=VCF'")
                        return False
                    state = "MetadataHeader"

                elif state == "MetadataHeader":

                    # Check format of metadata header section
                    if line.startswith("#CHROM"):
                        # Check to see that required fields are present in header line
                        line = line.lstrip("#")
                        line = line.split()

                        # Get the number of samples in the VCF
                        num_samples = len(line) - len(VCFHelper.HEADER_FIELDS)

                        expected_header_string  = "".join(VCFHelper.HEADER_FIELDS)
                        actual_header_string    = "".join(line[0:len(VCFHelper.HEADER_FIELDS)])

                        # Check to make sure the header fields are present in the correct order
                        if actual_header_string != expected_header_string:
                            logging.error("Invalid VCF file format! Invalid static column lables.\n"
                                          "Expected: %s\n"
                                          "Received: %s" % (expected_header_string, actual_header_string))
                            return False

                        # Update current state
                        state = "FirstRecord"

                    elif not line.startswith("##"):
                        # Return false if header line doesn't start with '##'
                        logging.error("Invalid VCF file format! Metadata header lines must begin with '##'!")
                        return False

                elif state == "FirstRecord":

                    # Validate first record
                    line = line.split()
                    if len(line) != num_samples + len(VCFHelper.HEADER_FIELDS):
                        logging.error("Invalid VCF file format! First record did not contain the correct number of columns!")
                        return False

                    # Return True as VCF has been validated to the point that I give a shit about
                    logging.debug("VCF file is valid: %s" % path)
                    return True

    @staticmethod
    def get_annotation_type(path):
        # Determine whether VCF file has been annotated by snpeff, annovar, or other
        with open(path, 'r') as vcf_fh:
            vcf_parser = vcf.Reader(vcf_fh)
            if "ANNOVAR_DATE" in vcf_parser.infos.keys():
                # Return 'annovar' if 'ANNOVAR_DATE' appears in VCF
                return VCFAnnotationType.ANNOVAR

            elif "ANN" in vcf_parser.infos.keys():
                # See if 'ANN' field with description 'Functional Annotations' is present
                desc = vcf_parser.infos["ANN"].desc
                if desc.startswith("Functional annotations"):
                    return VCFAnnotationType.SNPEFF
            else:
                return VCFAnnotationType.UNKNOWN

    @staticmethod
    def get_info_field_names(vcf_parser):
        # Return all available annotation fields for a VCF file (in order)
        # Get available field names
        field_names = OrderedDict.fromkeys(vcf_parser.infos.keys()).keys()
        if "ANN" in vcf_parser.infos.keys():
            # SNPEff annotated VCF file
            snpeff_fields = vcf_parser.infos["ANN"].desc.split("'")[1].split("|")
            snpeff_fields = [x.strip() for x in snpeff_fields]

            # Remove 'ANN' and replace with unpacked values
            ann_index = field_names.index("ANN")
            field_names = field_names[0:ann_index] + snpeff_fields + field_names[ann_index+1:]

        return field_names

    @staticmethod
    def get_snpeff_ann_field_names(vcf_parser):
        # Return all available annotation fields contained in the SNPeff 'ANN' section
        snpeff_fields = vcf_parser.infos["ANN"].desc.split("'")[1].split("|")
        return [x.strip() for x in snpeff_fields]
