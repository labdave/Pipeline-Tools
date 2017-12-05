import logging
import vcf
from collections import OrderedDict

class VCFRecoder:

    # Class for parsing VCF objects
    def __init__(self, vcf_file, min_call_depth=10, missing_data_char='.'):

        # Path to vcf file
        self.vcf_file = vcf_file

        # Mininum number of reads required to make a variant call
        self.min_call_depth = min_call_depth

        # Placeholder to use when data is missing
        self.missing_data_char = missing_data_char

    def recode_vcf(self):
        # Parse VCF and recode genotypes and output information as tab-delimited file

        logging.info("Reading vcf file: %s" % self.vcf_file)
        vcf_parser = vcf.Reader(open(self.vcf_file, "r"))

        # Get header columns
        header_data = self.get_header_data(vcf_parser)

        # Print header to stdout
        print "\t".join(header_data["header"])

        # Iterate over VCF records
        for record in vcf_parser:

            # Get basic information strings
            chrom   = record.CHROM
            pos     = record.POS
            db_id   = record.ID if record.ID is not None else self.missing_data_char
            ref     = record.REF

            # Get string representation of alternate allele
            alt     = self.get_alt_allele(record)
            qual    = record.QUAL
            filter  = record.FILTER if record.filter is not None else "."

            # Get ordered dictionary of variant info
            info        = get_variant_info(record)

            # Get ordered dictionary of sample genotypes
            genotypes   = get_genotypes(record.samples)

            # Iterate over samples genotypes at this locus
            for sample in record.samples:

                # Recode the genotype
                print sample

            # Print out a tab-delimited string version of the variant with the new genotypes
            print record

    def get_header_data(self, vcf_parser):
        # Return list of columns that will appear in the recoded VCF

        header_data = dict()

        # Add basic informational columns
        header_data["fixed"] = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"]

        # Add info columns parsed from metadata
        header_data["info"] = OrderedDict.fromkeys(vcf_parser.infos.keys()).keys()
        header_data["info"].remove("ALLELE_END")

        # Add sample names
        header_data["sample"] = vcf_parser.samples

        # Create entire header
        header_data["header"] = header_data["fixed"] + header_data["info"] + header_data["sample"]

        return header_data

    def recode_genotype(self, sample):
        pass

    def stringify_vcf_record(self, record):
        pass

    def get_alt_allele(self, record):
        # Return string representation of alternate allele

        # Case: One alternate allele
        if len(record.ALT) == 1 and record.ALT[0] is not None:
            return record.ALT[0]

        # Case: No alternate allele
        elif len(record.ALT) == 1 and record.ALT[0] is None:
            return self.missing_data_char

        elif len(record.ALT) > 1:
            # Return multiple alternate alleles as comma-separated list
            logging.warning("(VCFRecoder) Multiple alternate alleles detected for record: %s" % record)
            return ",".join([str(x) for x in record.ALT])

    def get_variant_info(self, record):
        # Return a



