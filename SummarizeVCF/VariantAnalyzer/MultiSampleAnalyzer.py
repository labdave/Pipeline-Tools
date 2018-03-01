from VCF import VCFHelper, VCFAnnotationType

class MultiSampleAnalyzer(object):

    def __init__(self, annotation_parser):
        # Class for parsing information from annotation fields in a standardized way
        self.annotation_parser = annotation_parser

    def declare_required_count_fields(self):

        required_fields = ["Missing GT", "Called GT", "Variant GT", "Heterozygous", "Homozygous-Alt", "Deletions", "Insertions",
                           "Monomorphs", "SNPs", "Ts", "Tv"]

        # Add nucleotide transitions
        for ref_base in ["A", "C", "G", "T"]:
            for alt_base in ["A", "C", "G", "T"]:
                if ref_base != alt_base:
                    required_fields.append("%s%s" % (ref_base, alt_base))

        # Add other stuff
        required_fields.append("dbSNP")

        avail_info = self.annotation_parser.get_available_info_fields()
        annotation_type = self.annotation_parser.get_annotation_type()

        if "Annotation_Impact" in avail_info:
            # Add SnpEff variant impact counts if they're present
            required_fields += ["%s_Impact" % x for x in ["MODIFIER", "LOW", "MODERATE", "HIGH"]]

        if annotation_type is VCFAnnotationType.SNPEFF:
            # Add notable SnpEff mutation types SnpEff
            required_fields += ["intergenic_region", "intron_variant",
                                "synonymous_variant","missense_variant",
                                "stop_gained","stop_lost",
                                "inframe_insertion","inframe_deletion",
                                "frameshift_variant"]

        elif annotation_type is VCFAnnotationType.ANNOVAR:
            required_fields += ["intergenic","intronic",
                                "synonymous_SNV", "nonsynonymous_SNV",
                                "stoploss","stopgain",
                                "nonframeshift_deletion","nonframeshift_insertion",
                                "frameshift_deletion","frameshift_insertion"]

        return required_fields

    def analyze(self, record, vcf_summary):
        # Process VCF record and add any necessary information to the VCFSummary

        # Determine whether variant is a spanning deletion (alternate allele will be '*')
        spanning_deletion = record.alleles[1] == "*"

        # Parse available variant annotations
        record_info = self.annotation_parser.get_info(record)

        # Get variant information before processing genotypes
        if record.is_indel or spanning_deletion:
            # Determine deletion length
            indel_len = abs(VCFHelper.get_variant_size(record))

        elif record.is_snp:
            # Determine the type of transition (e.g. A->G, G-T)
            transition_type = VCFHelper.get_snp_transition_type(record)
            ts_tv_status = "Ts" if record.is_transition else "Tv"

        variant_class   = VCFHelper.get_variant_class(record_info)
        is_dbsnp        = VCFHelper.is_dbsnp(record, record_info)
        variant_effect  = VCFHelper.get_snpeff_impact(record_info)

        # Process each sample
        for sample in record.samples:

            sample_name = sample.sample

            # Get sample genotype
            gt = record.genotype(sample_name)
            gt_type = gt.gt_type

            if gt_type is None:
                # Increment missing genotype count
                vcf_summary.add_count(sample_name, "Missing GT")
                continue

            # Increment called genotype count regardless of whether it's called a variant
            vcf_summary.add_count(sample_name, "Called GT")

            # Report sequencing depth for all called genotypes
            if hasattr(gt.data,"AD"):
                vcf_summary.add_depth(sample_name, depth=sum(gt.data.AD))

            # Report genotype quality for all called genotypes
            if hasattr(gt.data, "GQ"):
                vcf_summary.add_qual(sample_name, qual=gt.data.GQ)

            if not gt_type > 0:
                # Skip homo ref genotypes
                continue

            # Increment total number of variants for sample
            vcf_summary.add_count(sample_name, "Variant GT")

            # Report alternate allele frequency for all variants
            vcf_summary.add_aaf(sample_name, record.aaf[0])

            # Add information for whether variant is hetero or homo alternate allele
            if gt_type == 1:
                # Heterozygous
                vcf_summary.add_count(sample_name, "Heterozygous")

            else:
                # Homozygous alternate
                vcf_summary.add_count(sample_name, "Homozygous-Alt")

            # Add information for type of variant
            if record.is_indel and record.is_deletion or spanning_deletion:
                # Record deletion
                vcf_summary.add_count(sample_name, "Deletions")
                vcf_summary.add_del(sample_name, del_len=indel_len)

            elif record.is_indel:
                # Record deletion
                vcf_summary.add_count(sample_name, "Insertions")
                vcf_summary.add_ins(sample_name, ins_len=indel_len)

            elif record.is_snp:
                # Record SNP
                vcf_summary.add_count(sample_name, "SNPs")
                vcf_summary.add_count(sample_name, transition_type)
                vcf_summary.add_count(sample_name, ts_tv_status)

            elif record.is_monomorphic:
                # Record monomorphic loci
                vcf_summary.add_count(sample_name, "Monomorphs")

            elif record.is_sv:
                vcf_summary.add_count(sample_name, "Structural Variants")

            else:
                # Record unkown variant type
                vcf_summary.add_count(sample_name, "Unknown Variant Type")


            if is_dbsnp:
                # Increment if variant annotated as appearing in dbSNP
                vcf_summary.add_count(sample_name, "dbSNP")

            # Add additional information
            if variant_class is not None:
                vcf_summary.add_count(sample_name, variant_class)

            # Add snpeff variant effect info
            if variant_effect is not None:
                vcf_summary.add_count(sample_name, "%s_Impact" % variant_effect)
