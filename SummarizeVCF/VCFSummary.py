from SummarizeVCF.SampleSummary import SampleSummary

class VCFSummary(object):
    def __init__(self, sample_names, required_count_names=None, **kwargs):

        # List of samples to be contained in report
        self.sample_names   = sample_names

        # List of colnames that will be counted in report
        self.count_names    = [] if required_count_names is None else required_count_names

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

    def init_count_field(self, count_name):
        # Initialize a new count field for all samples
        for sample_name in self.sample_names:
            self.summary[sample_name].init_count(count_name)
        self.count_names.append(count_name)

    def add_count(self, sample_name, count_name):
        if count_name not in self.count_names:
            self.init_count_field(count_name)
        self.summary[sample_name].add_count(count_name)

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
        return self.__get_count_data_string() + self.__get_depth_summary_string() + self.__get_qual_summary_string() + self.__get_afs_summary_string() + self.__get_insert_summary_string() + self.__get_delete_summary_string()

    def __get_count_data_string(self):
        to_ret = "#COUNTS\n"
        count_names = [x.replace(" ", "_") for x in self.count_names]
        to_ret += "\t".join(["Sample"] + count_names) + "\n"
        for sample in self.sample_names:
            sample_data = [sample] + ["%d" % self.summary[sample].get_count(count_name) for count_name in self.count_names]
            to_ret += "\t".join(sample_data) + "\n"
        return to_ret

    def __get_depth_summary_string(self):
        to_ret = "#DEPTH\n"
        first_sample = True
        for sample in self.sample_names:
            labels, data = self.summary[sample].get_depth_summary()
            if first_sample:
                to_ret += "\t".join(["Sample"] + labels) + "\n"
                first_sample = False
            sample_data = [sample] + ["%d" % x for x in data]
            to_ret += "\t".join(sample_data) + "\n"
        return to_ret

    def __get_qual_summary_string(self):
        to_ret = "#QUAL\n"
        first_sample = True
        for sample in self.sample_names:
            labels, data = self.summary[sample].get_qual_summary()
            if first_sample:
                to_ret += "\t".join(["Sample"] + labels) + "\n"
                first_sample = False
            sample_data = [sample] + ["%d" % x for x in data]
            to_ret += "\t".join(sample_data) + "\n"
        return to_ret

    def __get_insert_summary_string(self):
        to_ret = "#INSERT_LEN\n"
        first_sample = True
        for sample in self.sample_names:
            labels, data = self.summary[sample].get_ins_summary()
            if first_sample:
                to_ret += "\t".join(["Sample"] + labels) + "\n"
                first_sample = False
            sample_data = [sample] + ["%d" % x for x in data]
            to_ret += "\t".join(sample_data) + "\n"
        return to_ret

    def __get_delete_summary_string(self):
        to_ret = "#DELETE_LEN\n"
        first_sample = True
        for sample in self.sample_names:
            labels, data = self.summary[sample].get_del_summary()
            if first_sample:
                to_ret += "\t".join(["Sample"] + labels) + "\n"
                first_sample = False
            sample_data = [sample] + ["%d" % x for x in data]
            to_ret += "\t".join(sample_data) + "\n"
        return to_ret

    def __get_afs_summary_string(self):
        to_ret = "#AAFS\n"
        first_sample = True
        for sample in self.sample_names:
            labels, data = self.summary[sample].get_aaf_summary()
            if first_sample:
                to_ret += "\t".join(["Sample"] + labels) + "\n"
                first_sample = False
            sample_data = [sample] + ["%d" % x for x in data]
            to_ret += "\t".join(sample_data) + "\n"
        return to_ret
