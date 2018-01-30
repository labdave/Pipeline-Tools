class SampleSummary(object):
    # Class for holding summary data about a VCF file
    def __init__(self, max_depth=500, max_qual=500, max_indel_len=250, num_afs_bins=20):
        self.summary        = {}
        self.max_depth      = max_depth
        self.max_qual       = max_qual
        self.max_indel_len  = max_indel_len
        self.num_afs_bins   = num_afs_bins
        self.depths         = [0] * (max_depth + 1)
        self.quals          = [0] * (max_qual + 1)
        self.inserts        = [0] * (max_indel_len + 1)
        self.deletes        = [0] * (max_indel_len + 1)
        self.afs            = [0] * num_afs_bins

    def add_count(self, key):
        # Increment a count by one
        if key in self.summary:
            self.summary[key] += 1
        else:
            self.summary[key] = 1

    def add_depth(self, depth):
        # Add variant depth to depth counts
        if depth >= self.max_depth:
            self.depths[self.max_depth] += 1
        else:
            self.depths[depth] += 1

    def add_qual(self, qual):
        # Add variant quality to qual counts
        if qual >= self.max_qual:
            self.quals[self.max_qual] += 1
        else:
            self.quals[qual] += 1

    def add_del(self, del_len):
        # Add deletion length to deletion length counts
        if del_len >= self.max_indel_len:
            self.deletes[self.max_indel_len] += 1
        else:
            self.deletes[del_len] += 1

    def add_ins(self, ins_len):
        # Add insertion length to insertion length counts
        if ins_len >= self.max_indel_len:
            self.inserts[self.max_indel_len] += 1
        else:
            self.inserts[ins_len] += 1

    def add_aaf(self, aaf):
        # Add alternative allele frequency (AAF) to AAF counts bin
        # Determine bin membership and increment appropriate AFS count
        aaf_bin = int(aaf//self.num_afs_bins)
        self.afs[aaf_bin] += 1
