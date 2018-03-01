class VCFSummaryParser:

    def __init__(self, vcf_summary_file):
        self.vcf_summary_file = vcf_summary_file

    def parse(self):
        curr_section =