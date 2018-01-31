from SampleSummary import SampleSummary

class VCFSummary(object):
    def __init__(self, max_depth, max_qual, max_indel_len, num_afs_bins):

        # Upper bin for depth, qual, indel length counts
        self.max_depth      = max_depth
        self.max_qual       = max_qual
        self.max_indel_len  = max_indel_len

        # Number of bins to use for allele frequency spectrum counts
        self.num_afs_bins   = num_afs_bins

        # Main summary where values are incremented
        self.summary = {}

    def add_sample_summary(self, summary_name):
        self.summary[summary_name] = SampleSummary(self.max_depth, self.max_qual, self.max_indel_len, self.num_afs_bins)
