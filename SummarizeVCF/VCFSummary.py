from SummarizeVCF.SampleSummary import SampleSummary

class VCFSummary(object):
    def __init__(self, sample_names, count_names, **kwargs):

        # List of samples to be contained in report
        self.sample_names   = sample_names

        # List of colnames that will be counted in report
        self.count_names    = count_names

        # Upper bin for depth, qual, indel distribution summaries
        self.max_indel_len  = kwargs.get("max_indel_len",   100)
        self.max_depth      = kwargs.get("max_depth",       500)
        self.max_qual       = kwargs.get("max_qual",        250)

        # Number of bins to use for allele frequency spectrum counts
        self.num_afs_bins   = kwargs.get("num_afs_bins",    20)

        # Main summary where values are incremented
        self.summary = self.init_summary()

    def init_summary(self):
        # Create new summaries for each sample and initialize counters for each count_name
        summary = {}
        for sample_name in self.sample_names:
            summary[sample_name] = SampleSummary(self.max_depth, self.max_qual, self.max_indel_len, self.num_afs_bins)
            for count_name in self.count_names:
                summary[sample_name].init_count(count_name)
        return summary

    def add_count(self, sample_name, key):
        self.summary[sample_name].add_count(key)

    def add_depth(self, sample_name, depth):
        self.summary[sample_name].add_depth(depth)

    def add_qual(self, sample_name, qual):
        self.summary[sample_name].add_qual(qual)

    def add_del(self, sample_name, del_len):
        self.summary[sample_name].add_del(del_len)

    def add_ins(self, sample_name, ins_len):
        self.summary[sample_name].add_ins(ins_len)

    def add_aaf(self, sample_name, aaf):
        self.summary[sample_name].add_aaf(aaf)

    def __str__(self):
        # Get string representation of current report
        pass
