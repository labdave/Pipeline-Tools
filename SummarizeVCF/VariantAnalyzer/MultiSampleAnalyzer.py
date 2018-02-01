class MultiSampleAnalyzer(object):
    # Names of counters that need to be initialized in vcfsummary
    COUNTER_NAMES = []

    def __init__(self, annotation_parser):
        # Class for parsing information from annotation fields in a standardized way
        self.annotation_parser = annotation_parser

    @staticmethod
    def analyze(record, vcf_summary):
        # Process VCF record and add any necessary information to the VCFSummary

        # Get data members that would take a while to process
        if record.is_indel and record.is_deletion:
            # Determine deletion length
            indel_len = 0

        elif record.is_indel:
            # Determine insertion length
            indel_len = 0

        elif record.is_snp:
        # Determine the type of snp
            pass


        # Process each sample
        for sample in record.samples:

            sample_name = sample.name

            # Get sample genotype
            gt = record.genotype(sample_name)
            gt_type = gt.gt_type

            # Report sequencing depth for all samples, positions
            vcf_summary.add_depth(sample_name, depth=gt.data[2])

            if gt_type is None:
                # Increment missing genotype count
                vcf_summary.add_count(sample_name, "Missing GT")
                continue

            # Increment called genotype count regardless of whether it's called a variant
            vcf_summary.add_count(sample_name, "Called GT")

            if not gt_type > 0:
                # Skip homo ref genotypes
                continue

            # Increment total number of variants for sample
            vcf_summary.add_count(sample_name, "Variants")

            # Report variant call quality and alternate allele frequency
            vcf_summary.add_qual(sample_name, qual=gt.data[3])
            vcf_summary.add_aaf(sample_name, record.aaf)

            # Add information for whether variant is hetero or homo alternate allele
            if gt_type == 1:
                # Heterozygous
                vcf_summary.add_count(sample_name, "Heterozygotes")

            else:
                # Homozygous alternate
                vcf_summary.add_count(sample_name, "Homozygous-Alt")

            # Add information for type of variant
            if record.is_indel and record.is_deletion:
                # Record deletion
                vcf_summary.add_count(sample_name, "Deletions")
                vcf_summary.add_del(sample_name, del_len=indel_len)

            elif record.is_indel:
                # Record deletion
                vcf_summary.add_count(sample_name, "Insertions")
                length = 0
                vcf_summary.add_ins(sample_name, ins_len=indel_len)

            elif record.is_snp:
                # Record SNP
                vcf_summary.add_count(sample_name, "SNPs")

            elif record.is_monomorphic:
                # Record monomorphic loci
                vcf_summary.add_count(sample_name, "Monomorphic")

            else:
                # Record unkown variant type
                vcf_summary.add_count(sample_name, "Unknown Variant Type")