import logging
import vcf
from collections import OrderedDict

class VCFRecoder(object):

    # Class for parsing VCF objects
    def __init__(self, vcf_file, **kwargs):

        # Path to vcf file
        self.vcf_file = vcf_file

        # Mininum number of reads required to make a variant call
        self.min_call_depth = kwargs.get("min_call_depth", 10)

        # Placeholder to use when data is missing
        self.missing_data_char = kwargs.get("missing_data_char", '.')

        # Placeholder for uncalled genotypes
        self.missing_gt_char    = kwargs.get("missing_gt_char", 'NA')

        # Initialize class variables
        logging.debug("Initializing VCF parser for file: %s" % self.vcf_file)
        self.parser         = vcf.Reader(open(self.vcf_file, "r"))

        # Names of columns to include (in the order they will appear in the recoded VCF)
        self.fixed_columns  = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"]    # Fixed VCF columns
        self.info_columns   = self.get_info_field_names()                               # Analysis specific information columns
        self.sample_names   = self.parser.samples                                       # Names of samples included in VCF

    def recode_vcf(self):
        # Parse VCF and recode genotypes and output information as tab-delimited file

        # Combine columns into a single header
        colnames = self.fixed_columns + self.info_columns + self.sample_names
        num_cols = len(colnames)

        # Print header to stdout
        print "\t".join(colnames)

        # Iterate over VCF records
        for record in self.parser:

            # Check to make sure there is only one alternate allele for the record
            if len(record.ALT) > 1:
                # Raise error because it's just easier if we make people normalize their data beforehand
                logging.error("(VCFRecoder) Multiallelic error! All VCF records must have only one alternate allele. "
                              "Received the following record:\n%s" % record)
                # Failing to raise this error could cause some weirdo unintended bugs
                raise IOError("Multiple alternate alleles detected at one or more positions in vcf file!")

            record_data = self.get_fixed_data(record)

            record_data += self.get_info_data(record)

            record_data += self.get_genotype_data(record)

            if len(record_data) != num_cols:
                logging.error("(VCFRecoder) Record doesn't contain the same number of columns as header:\n%s" % record)
                raise IOError("VCF Record contains different number of columns than header!")

            # Combine the list into a string and print to stdout
            record_data = [str(x) for x in record_data]
            print "\t".join(record_data)

    def get_fixed_data(self, record):
        # Add fixed data columns first
        data = [record.CHROM, record.POS]

        # Add dbsnp id
        data.append(record.ID if record.ID is not None else self.missing_data_char)

        # Add information for reference base
        data.append(record.REF)

        # Add alt allele information
        data.append(self.missing_data_char if len(record.ALT) == 0 else record.ALT[0])

        # Add quality information
        data.append(record.QUAL)

        # Add filter information
        data.append(self.missing_data_char if record.FILTER is None else record.FILTER)
        return data

    def get_info_data(self, record):
        data = []
        # Loop through required info columns in order and get their value for this record
        for info_column in self.info_columns:
            # Make sure the info columns is actually contained in the record
            if info_column in record.INFO:
                info_data = record.INFO[info_column]

                if isinstance(info_data, list):
                    # Case: record.INFO[info_column] is a list of values

                    if info_data[0] is None and len(info_data) == 1:
                        # Case: Missing data
                        data.append(self.missing_data_char)
                    else:
                        # Case: One or more data points as a list
                        data.append(",".join([str(x) for x in info_data]))

                elif info_data is None:
                    data.append(self.missing_data_char)

                else:
                    # Case: Data is a single value
                    data.append(str(info_data))

            else:
                data.append(self.missing_data_char)
        return data

    def get_genotype_data(self, record):
        data = []
        for sample in self.sample_names:
            data.append(self.recode_genotype(record.genotype(sample)))
        return data

    def recode_genotype(self, sample_genotype):
        # Returns numerical representation of genotype call
        #  -1: Homozygous ref with >min_call_depth # reads supporting REF
        #   1: Variant called with >min_call_depth # reads supporting ALT
        #   0: <min_call_depth # reads supporting either REF or ALT for a sample

        if not sample_genotype.called:
            # Return 0 if sample not called
            return self.missing_gt_char

        elif sample_genotype.gt_type == 0:
            # Case: Called homozygous REF

            # Get allele depth of reference allele
            dp = sample_genotype.data.AD[0]

            if dp >= self.min_call_depth:
                # Case: DP > min_call_depth so report that we are confident GT is homo REF (-1)
                return "-1"
            else:
                # Case: DP < min_call_depth so report a decimal between 0 and -1 how close the supporting reads are to min_call_depth
                recoded = -1.0 * (dp/float(self.min_call_depth))
                if recoded == 0.0:
                    return "-0"
                return str(recoded)

        elif sample_genotype.gt_type != 0:
            # Case: Called hetero or homozygous ALT

            # Get allele depth of alternate allele
            dp = sample_genotype.data.AD[1]

            if dp >= self.min_call_depth:
                return "1"
            else:
                recoded = dp/float(self.min_call_depth)
                if recoded == 0.0:
                    return "0"
                else:
                    return str(recoded)

    def get_info_field_names(self):
        info_fields = OrderedDict.fromkeys(self.parser.infos.keys()).keys()
        logging.debug("(VCFRecoder) Getting following info fields:\n%s" % info_fields)
        return info_fields







