import logging

from VCFRecoder import VCFRecoder

class SnpEffVCFRecoder(VCFRecoder):

    def __init__(self, vcf_file, **kwargs):

        # Call super constructor
        super(SnpEffVCFRecoder, self).__init__(vcf_file, **kwargs)

    def get_info_field_names(self):
        # Return SnpEff annotation fields

        # Make sure SnpEff 'ANN' field detected in INFO section
        if not "ANN" in self.parser.infos:
            logging.error("Invalid SnpEff VCF! No 'ANN' annotation section declared in VCF INFO section header. "
                          "Are you sure this vcf was annotated with SnpEff?")
            raise IOError("Invalid SnpEff VCF!")

        return self.parser.infos["ANN"].desc.split("'")[1].split("|")

    def get_info_data(self, record):
        # Return only SnpEff annotations for a variant record
        data = record.INFO["ANN"][0].split("|")
        return [self.missing_data_char if x == "" else x for x in data]


