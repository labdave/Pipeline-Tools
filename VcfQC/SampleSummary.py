class SampleSummary(object):
    # Class for holding summary data about a VCF file
    def __init__(self, max_depth=500, max_qual=500, max_indel_len=250, num_afs_bins=20):

        # Main summary where values are incremented
        self.summary        = {}

        # Upper bin for depth, qual, indel length counts
        self.max_depth      = max_depth
        self.max_qual       = max_qual
        self.max_indel_len  = max_indel_len

        # Number of bins to use for allele frequency spectrum counts
        self.num_afs_bins   = num_afs_bins

        # Initialize count arrays for holding histogram counts for each data type
        self.depths         = [0] * (max_depth + 1)
        self.quals          = [0] * (max_qual + 1)
        self.inserts        = [0] * (max_indel_len + 1)
        self.deletes        = [0] * (max_indel_len + 1)
        self.afs            = [0] * num_afs_bins

    def add_count(self, key):
        # Increment a count in the main summary by one
        if key in self.summary:
            self.summary[key] += 1
        else:
            self.summary[key] = 1

    def add_depth(self, depth):
        # Increment depth histogram bin corresponding to observed depth
        if depth >= self.max_depth:
            self.depths[self.max_depth] += 1
        else:
            self.depths[depth] += 1

    def add_qual(self, qual):
        # Increment variant histogram bin corresponding to observed variant quality score
        if qual >= self.max_qual:
            self.quals[self.max_qual] += 1
        else:
            self.quals[qual] += 1

    def add_del(self, del_len):
        # Increment deletion length histogram bin corresponding to observed deletion
        if del_len >= self.max_indel_len:
            self.deletes[self.max_indel_len] += 1
        else:
            self.deletes[del_len] += 1

    def add_ins(self, ins_len):
        # Increment insertion length histogram bin corresponding to observed insertion
        if ins_len >= self.max_indel_len:
            self.inserts[self.max_indel_len] += 1
        else:
            self.inserts[ins_len] += 1

    def add_aaf(self, aaf):
        # Add alternative allele frequency (AAF) to AAF histogram
        # Determine bin membership
        aaf_bin = int(aaf//self.num_afs_bins)
        # Increment appropriate count
        self.afs[aaf_bin] += 1
