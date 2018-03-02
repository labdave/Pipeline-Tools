import logging

class VCFSummeryMergeError(Exception):
    pass

def merge_two_vcf_summaries(vcf_summary_1, vcf_summary_2):

    # Check that summaries can be merged
    can_merge, reason = can_merge_summaries(vcf_summary_1, vcf_summary_2)
    if not can_merge:
        logging.error("Unable to merge VCF Summaries! %s" % reason)
        raise VCFSummeryMergeError("Unable to merge VCF summaries!")

    # Init uninitialized counts in each VCF summary
    count_names_1 = vcf_summary_1.count_names
    count_names_2 = vcf_summary_2.count_names
    for count_name in count_names_1:
        # Add counts from summary 1 that don't appear in summary 2
        if count_name not in count_names_2:
            vcf_summary_2.init_count_field(count_name)
    for count_name in count_names_2:
        # Add counts from summary 2 that don't appear in summary 1
        if count_name not in count_names_1:
            vcf_summary_1.init_count_field(count_name)

    # Merge all sample summaries that appear in both VCFSummaries
    sample_names_1 = vcf_summary_1.sample_names
    sample_names_2 = vcf_summary_2.sample_names
    for sample_name in sample_names_1:
        if sample_name in sample_names_2:
            merge_sample_summaries(vcf_summary_1.summary[sample_name], vcf_summary_2.summary[sample_name])

    # Add any samples in VCFSummary_2 that don't appear in VCFSummary 1
    for sample_name in sample_names_2:
        if sample_name not in sample_names_1:
            vcf_summary_1.summary[sample_name] = vcf_summary_2.summary[sample_name]
            vcf_summary_1.sample_names.append(sample_name)

    # Return merged VCFSummary
    return vcf_summary_1

def merge_sample_summaries(sample_1, sample_2):

    # Combine counts
    for count_name in sample_1.summary:
        sample_1.summary[count_name] += sample_2.summary[count_name]

    # Combine depths
    sample_1.depths = [sample_1.depths[i] + sample_2.depths[i] for i in range(len(sample_1.depths))]

    # Combine depths
    sample_1.quals = [sample_1.quals[i] + sample_2.quals[i] for i in range(len(sample_1.quals))]

    # Combine depths
    sample_1.afs = [sample_1.afs[i] + sample_2.afs[i] for i in range(len(sample_1.afs))]

    # Combine depths
    sample_1.inserts = [sample_1.inserts[i] + sample_2.inserts[i] for i in range(len(sample_1.inserts))]

    # Combine depths
    sample_1.deletes = [sample_1.deletes[i] + sample_2.deletes[i] for i in range(len(sample_1.deletes))]

def can_merge_summaries(vcf_summary_1, vcf_summary_2):
    # Check to see if two summaries can be merged
    can_merge = True
    reason = ""
    if vcf_summary_1.max_indel_len != vcf_summary_2.max_indel_len:
        # Check max indel length the same
        can_merge = False
        reason = "Differing max indel lengths: %s vs. %s!" % (vcf_summary_1.max_indel_len, vcf_summary_2.max_indel_len)
    elif vcf_summary_1.max_depth != vcf_summary_2.max_depth:
        # Check max depth summary the same
        can_merge = False
        reason = "Differing max depths: %s vs. %s!" % (vcf_summary_1.max_depth, vcf_summary_2.max_depth)
    elif vcf_summary_1.max_qual != vcf_summary_2.max_qual:
        # Check max quality summary the same
        can_merge = False
        reason = "Differing max qual: %s vs. %s!" % (vcf_summary_1.max_qual, vcf_summary_2.max_qual)
    elif vcf_summary_1.num_afs_bins != vcf_summary_2.num_afs_bins:
        # Check same number of AAFS bins
        can_merge = False
        reason = "Differing num AAFS bins: %s vs. %s!" % (vcf_summary_1.num_afs_bins, vcf_summary_2.num_afs_bins)
    return can_merge, reason


