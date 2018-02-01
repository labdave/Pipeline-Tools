import logging
import vcf

from VCF import AnnotationParser
from SummarizeVCF.VariantAnalyzer import VariantAnalyzerFactory
from VCFSummary import VCFSummary

class VCFSummarizer(object):

    # Class for parsing VCF objects
    def __init__(self, vcf_file, summary_type, **kwargs):

        # Get VCFParser
        self.vcf_parser = vcf.Reader(open(vcf_file, "r"))

        # Get annotation parser specific to this vcf
        self.annotation_parser = AnnotationParser(self.vcf_parser)

        # Get variant analyzer
        variant_analyzer_factory = VariantAnalyzerFactory(self.annotation_parser)
        self.variant_analyzer = variant_analyzer_factory.get_analyzer(summary_type)

        # Create VCFSummary
        self.summary = VCFSummary(sample_names=self.vcf_parser.samples,
                                  count_names=self.variant_analyzer.COUNTER_NAMES,
                                  **kwargs)

    def summarize(self):
        # Parse VCF and generate summary statistics

        # Num lines processed
        processed = 0

        for variant_record in self.vcf_parser:

            if processed % 100000 == 0:
                logging.info("(VCFSummarizer) Processed %d records!" % processed)

            # Check to make sure there is only one alternate allele for the record
            if len(variant_record.ALT) > 1:
                # Raise error because it's just easier if we make people normalize their data beforehand
                logging.error("(VCFSummarizer) Multiallelic error! All VCF records must have only one alternate allele. "
                              "Received the following record:\n%s" % variant_record)

                # Failing to raise this error could cause some weirdo unintended bugs
                raise IOError("Multiple alternate alleles detected at one or more positions in vcf file!")

            # Add data from summarizer to report
            self.variant_analyzer.analyze(record=variant_record, vcf_summary=self.summary)










            #print record.aaf
            #print record.get_hets()
            #print len(record.get_unknowns())
            #print record.var_subtype

            # Increment number of records processed
            processed += 1

    def get_summary(self):
        return self.summary







