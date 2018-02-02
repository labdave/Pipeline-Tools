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

    @staticmethod
    def get_variant_size(record):
        # Return length of an insertion or deletion (Assumes two alleles)
        return  len(record.alleles[1]) - len(record.alleles[0])

    @staticmethod
    def get_snp_transition_type(record):
        return "%s%s" % (record.alleles[0], record.alleles[1])

    @staticmethod
    def get_variant_class(record_info):
        # Return the type of variant (e.g. intronic, missense)

        if "ExonicFunc.refGene" in record_info and record_info["ExonicFunc.refGene"] is not None:
            # See if annovar annotation is present for coding variant
            # Remove any commas
            return record_info["ExonicFunc.refGene"].split("\\x3")[0]

        elif "Func.refGene" in record_info and record_info["Func.refGene"] is not None:
            # Check to see if annovar annotatio is present for intergenic/intronic variants
            return record_info["Func.refGene"].split("\\x3")[0]

        elif "Annotation" in record_info and record_info["Annotation"] is not None:
            # See if SNPeff annotation is present
            # Return only first annotation (multiple annotations can be joined by '&')
            return record_info["Annotation"].split("&")[0]

        return None

    @staticmethod
    def get_snpeff_impact(record_info):
        # Return the predicted effect of the mutation
        if "Annotation_Impact" in record_info:
            return record_info["Annotation_Impact"]
        return None

    @staticmethod
    def is_dbsnp(record, record_info):

        # Check to see if rsID exists
        if record.ID is not None and record.ID.startswith("rs"):
            return True

        # Check to see if Annovar annotations have dbSNP id
        for info in record_info:
            if info.startswith("snp1") and record_info[info] is not None:
                if record_info[info].startswith("rs"):
                    return True

        return False

