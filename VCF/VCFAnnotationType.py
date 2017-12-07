class VCFAnnotationType:
    # Static class for holding reserved VCF annotation types
    ANNOVAR = "ANNOVAR"
    SNPEFF  = "SNPEFF"
    UNKNOWN = "OTHER"
    annotation_types = [ANNOVAR, SNPEFF, UNKNOWN]
