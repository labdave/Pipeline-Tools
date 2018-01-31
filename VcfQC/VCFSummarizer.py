import logging
from collections import OrderedDict
import vcf

from SampleSummary import SampleSummary
from VCFSummary import VCFSummary

class VCFSummarizer(object):

    # Class for parsing VCF objects
    def __init__(self, vcf_file, out_file, **kwargs):

        # Path to vcf file
        self.vcf_file = vcf_file

        # Path to recoded output file
        self.out_file = out_file

        # Missing data char
        self.missing_data_char = kwargs.get("missing_data_char", ".")

        # Initialize VCF parser
        logging.debug("Initializing VCF parser for file: %s" % self.vcf_file)
        self.parser = vcf.Reader(open(self.vcf_file, "r"))

        # Get list of info columns
        self.info_columns = self.get_info_field_names()
        logging.debug("INFO columns:\n%s" % ", ".join(self.info_columns))

        # Get sample names
        self.sample_names   = self.parser.samples
        logging.debug("Samples:\n%s" % ", ".join(self.sample_names))

        # Initialize VCF summary object
        max_indel_len   = kwargs.get("max_indel_len", 100)
        max_depth       = kwargs.get("max_depth", 500)
        max_qual        = kwargs.get("max_qual", 250)
        num_afs_bins    = kwargs.get("num_afs_bins", 20)
        self.summary    = VCFSummary(max_indel_len, max_depth, max_qual, num_afs_bins)

    def summarize_vcf(self):
        # Parse VCF and generate summary statistics

        # Num lines processed
        processed = 0

        for record in self.parser:

            if processed % 100000 == 0:
                logging.info("(VCFSummarizer) Processed %d records!" % processed)

            # Check to make sure there is only one alternate allele for the record
            if len(record.ALT) > 1:
                # Raise error because it's just easier if we make people normalize their data beforehand
                logging.error("(VCFSummarizer) Multiallelic error! All VCF records must have only one alternate allele. "
                              "Received the following record:\n%s" % record)
                # Failing to raise this error could cause some weirdo unintended bugs
                raise IOError("Multiple alternate alleles detected at one or more positions in vcf file!")

            # Process missing genotype
            for unknown in record.get_unknowns():
                self.process_unknown_genotype(unknown)

            # Process called heterozygotes
            for het in record.get_hets():
                self.process_het_genotype(het)
                #self.summary[het.sample].add_count("heterozygotes")
                #self.summary[het.sample].add_depth(het[1])

            # Process called homozygous alternate genotypes
            for hom_alt in record.get_hom_alts():
                self.process_hom_alt_genotype(hom_alt)

            # Process called homozygous reference genotypes
            for hom_ref in record.get_hom_refs():
                self.process_hom_ref_genotype(hom_ref)

            print record.aaf
            print record.get_hets()
            print len(record.get_unknowns())
            print record.var_subtype

            # Increment number of records processed
            processed += 1

    def process_unknown_genotype(self, call_obj):
        pass

    def process_het_genotype(self, genotype):
        pass

    def process_hom_alt_genotype(self, genotype):
        pass

    def process_hom_ref_genotype(self, genotype):
        pass

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
            data.append(self.get_recoded_genotype(record.genotype(sample)))
        return data

    def get_recoded_genotype(self, sample_genotype):
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







