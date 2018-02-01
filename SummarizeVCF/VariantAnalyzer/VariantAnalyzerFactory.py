import logging

from MultiSampleAnalyzer import MultiSampleAnalyzer

class VariantAnalyzerFactorError(Exception):
    pass

class VariantAnalyzerFactory:
    # Factory class for initializing variant analysis types

    # Available summary types
    SUMMARY_TYPES = {"Multisample" : MultiSampleAnalyzer}

    def __init__(self, annotation_parser):
        self.annotation_parser = annotation_parser

    def get_analyzer(self, summary_type):
        # Returns a VariantAnalyzer instance created with the existing annotation parser
        if summary_type not in VariantAnalyzerFactory.SUMMARY_TYPES:
            logging.error("(VariantAnalyzerFactory) Cannot create variant analyzer for summary type '%s'" % summary_type)
            logging.error("(VariantAnalyzerFactory) Supported types: %s" % ", ".join(VariantAnalyzerFactory.SUMMARY_TYPES.keys()))
            raise VariantAnalyzerFactorError("Invalid VariantAnalyzer type!")
        return self.SUMMARY_TYPES[summary_type](self.annotation_parser)

    @staticmethod
    def get_available_types():
        return VariantAnalyzerFactory.SUMMARY_TYPES