from VcfQC.SampleSummary import SampleSummary

class VCFSummary(object):
    def __init__(self, sample_names, max_depth, max_qual, max_indel_len, num_afs_bins):

        # Names of samples in VCFSummary
        self.sample_names = sample_names

        # Upper bin for depth, qual, indel length counts
        self.max_depth      = max_depth
        self.max_qual       = max_qual
        self.max_indel_len  = max_indel_len

        # Number of bins to use for allele frequency spectrum counts
        self.num_afs_bins   = num_afs_bins

        # Main summary where values are incremented
        self.summary = self.init_summary()

    def init_summary(self):
        summary = {}
        for sample_name in self.sample_names:
            summary[sample_name] = SampleSummary(self.max_depth, self.max_qual, self.max_indel_len, self.num_afs_bins)
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

    def process_record(self, record):

        # Get data members that would take a while to process
        if record.is_indel and record.is_deletion:
            # Determine deletion length
            indel_len = 0

        elif record.is_indel:
            # Determine insertion length
            indel_len = 0

        elif record.is_snp:
        # Determine the type of snp


        # Process each sample
        for sample_name in self.sample_names:

            # Get sample genotype
            gt = record.genotype(sample_name)
            gt_type = gt.gt_type

            # Report sequencing depth for all samples, positions
            self.summary.add_depth(sample_name, depth=gt.data[2])

            if gt_type is None:
                # Increment missing genotype count
                self.summary.add_count(sample_name, "Missing GT")
                continue

            # Increment called genotype count regardless of whether it's called a variant
            self.summary.add_count(sample_name, "Called GT")

            if not gt_type > 0:
                # Skip homo ref genotypes
                continue

            # Increment total number of variants for sample
            self.summary.add_count(sample_name, "Variants")

            # Report variant call quality and alternate allele frequency
            self.summary.add_qual(sample_name, qual=gt.data[3])
            self.summary.add_aaf(sample_name, record.aaf)

            # Add information for whether variant is hetero or homo alternate allele
            if gt_type == 1:
                # Heterozygous
                self.summary.add_count(sample_name, "Heterozygotes")

            else:
                # Homozygous alternate
                self.summary.add_count(sample_name, "Homozygous-Alt")

            # Add information for type of variant
            if record.is_indel and record.is_deletion:
                # Record deletion
                self.summary.add_count(sample_name, "Deletions")
                self.summary.add_del(sample_name, del_len=indel_len)

            elif record.is_indel:
                # Record deletion
                self.summary.add_count(sample_name, "Insertions")
                length = 0
                self.summary.add_ins(sample_name, ins_len=indel_len)

            elif record.is_snp:
                # Record SNP
                self.summary.add_count(sample_name, "SNPs")

            elif record.is_monomorphic:
                # Record monomorphic loci
                self.summary.add_count(sample_name, "Monomorphic")

            else:
                # Record unkown variant type
                self.summary.add_count(sample_name, "Unknown Variant Type")
