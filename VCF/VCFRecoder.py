import logging
from collections import OrderedDict
import vcf

class VCFRecoder(object):

    # Class for parsing VCF objects
    def __init__(self, vcf_file, out_file, **kwargs):

        # Path to vcf file
        self.vcf_file = vcf_file

        # Path to recoded output file
        self.out_file = out_file

        # Mininum number of reads required to make a variant call
        self.min_call_depth             = kwargs.get("min_call_depth",      10)

        # Placeholder to use when data is missing
        self.missing_data_char          = kwargs.get("missing_data_char",   '.')

        # Placeholder for uncalled genotypes
        self.missing_gt_char            = kwargs.get("missing_gt_char",     'NA')

        # Specify whether or not to allow more than one multiple allele
        self.multiallelic               = kwargs.get("multiallelic",        False)

        # Path to file specifying info columns to include
        # If not specified, all info columns will be included
        self.info_column_file           = kwargs.get("info_column_file",        None)

        # Initialize class variables
        logging.debug("Initializing VCF parser for file: %s" % self.vcf_file)
        self.parser         = vcf.Reader(open(self.vcf_file, "r"))

        # Names of columns to include (in the order they will appear in the recoded VCF)
        self.fixed_columns  = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"]    # Names of columns that are fixed across VCF files
        self.info_columns   = self.get_info_field_names()                               # Names of INFO columns declared in metadata section
        self.sample_names   = self.parser.samples                                       # Names of samples included in current VCF

        # Parse file to get names of INFO columns to include in output
        self.info_to_include    = self.__parse_info_to_include_file()

        # Indices of INFO columns to include in output
        self.info_include_cols  = None

        # Subset info to include if selecting only certain INFO columns
        logging.debug("INFO columns available in VCF:\n%s" % self.info_columns)
        if self.info_to_include is not None:
            # Check to make sure columns requested in info-column-file actually exist in VCF
            self.__check_info_to_include()
            self.info_include_cols = [self.info_columns.index(i) for i in self.info_to_include]

            # Report columns actually being output to recoded VCF file
            logging.debug("INFO columns to be included in output: %s" % [self.info_columns[i] for i in self.info_include_cols])

    def recode_vcf(self):
        # Parse VCF and recode genotypes and output information as tab-delimited file

        # Open output file for writing
        out_file_handle = open(self.out_file, "w")

        # Combine columns into a single header
        if self.info_to_include is None:
            # Case: Get all INFO columns
            colnames = self.fixed_columns + self.info_columns + self.sample_names
        else:
            # Case: Include only certain INFO columns
            colnames = self.fixed_columns + [self.info_columns[i] for i in self.info_include_cols] + self.sample_names

        # Total number of data columns to get for each VCF record
        num_cols = len(colnames)

        # Write header to file
        out_file_handle.write("%s\n" % "\t".join(colnames))

        # Iterate over VCF records
        for record in self.parser:

            # Check to make sure there is only one alternate allele for the record
            if len(record.ALT) > 1 and not self.multiallelic:
                # Raise error because it's just easier if we make people normalize their data beforehand
                logging.error("(VCFRecoder) Multiallelic error! All VCF records must have only one alternate allele. "
                              "Received the following record:\n%s" % record)
                # Failing to raise this error could cause some weirdo unintended bugs
                raise IOError("Multiple alternate alleles detected at one or more positions in vcf file!")

            # Get fixed data
            record_data = self.get_fixed_data(record)

            # Get info column data
            if self.info_to_include is None:
                # Don't subset info columns
                record_data += self.get_info_data(record)
            else:
                # Subset info columns
                record_data += [self.get_info_data(record)[i] for i in self.info_include_cols]

            # Add recoded genotypes for each sample
            record_data += self.get_genotype_data(record)

            if len(record_data) != num_cols:
                logging.error("(VCFRecoder) Record doesn't contain the same number of columns as header:\n%s" % record)
                raise IOError("VCF Record contains different number of columns than header!")

            # Combine the list into a string and print to stdout
            record_data = [str(x) for x in record_data]

            # Write to file
            out_file_handle.write("%s\n" % "\t".join(record_data))

        # Close output file
        out_file_handle.close()

    def get_fixed_data(self, record):
        # Add fixed data columns first
        data = [record.CHROM, record.POS]

        # Add dbsnp id
        data.append(record.ID if record.ID is not None else self.missing_data_char)

        # Add information for reference base
        data.append(record.REF)

        # Add alt allele information
        data.append(self.missing_data_char if len(record.ALT) == 0 else record.ALT[0])

        # Add quality information
        data.append(record.QUAL if record.QUAL is not None else self.missing_data_char)

        # Add filter information
        filter = self.missing_data_char if record.FILTER is None else record.FILTER

        # Do additional check (Really only used for MUTECT output)
        if isinstance(filter, list):
            if len(filter) == 0:
                filter = "PASSED"
            else:
                filter = ",".join(filter)
        data.append(filter)

        # Return fixed data values
        return data

    def get_info_data(self, record):
        data = []
        # Loop through required info columns in order and get their value for this record
        for info_column in self.info_columns:
            # Make sure the info columns is actually contained in the record
            if info_column in record.INFO:
                info_data = record.INFO[info_column]

                if isinstance(info_data, list):
                    # Case: record.INFO[info_column] is a list of values

                    if info_data[0] is None and len(info_data) == 1:
                        # Case: Missing data
                        data.append(self.missing_data_char)
                    else:
                        # Case: One or more data points as a list
                        data.append(",".join([str(x) for x in info_data]))

                elif info_data is None:
                    data.append(self.missing_data_char)

                else:
                    # Case: Data is a single value
                    data.append(str(info_data))

            else:
                data.append(self.missing_data_char)
        return data

    def get_genotype_data(self, record):
        data = []
        for sample in self.sample_names:
            data.append(self.get_recoded_genotype(record.genotype(sample)))
        return data

    def get_recoded_genotype(self, sample_genotype):
        # Returns numerical representation of genotype call
        #  -1: Homozygous ref with >min_call_depth # reads supporting REF
        #   1: Variant called with >min_call_depth # reads supporting ALT
        #   0: <min_call_depth # reads supporting either REF or ALT for a sample

        if not sample_genotype.called:
            # Return 0 if sample not called
            return self.missing_gt_char

        elif sample_genotype.gt_type == 0:
            # Case: Called homozygous REF

            # Get allele depth of reference allele
            dp = sample_genotype.data.AD[0]

            if dp >= self.min_call_depth:
                # Case: DP > min_call_depth so report that we are confident GT is homo REF (-1)
                return "-1"
            else:
                # Case: DP < min_call_depth so report a decimal between 0 and -1 how close the supporting reads are to min_call_depth
                recoded = -1.0 * (dp/float(self.min_call_depth))
                if recoded == 0.0:
                    return "-0"
                return str(recoded)

        elif sample_genotype.gt_type != 0:
            # Case: Called hetero or homozygous ALT

            # Get allele depth of alternate allele
            dp = sample_genotype.data.AD[1]

            if dp >= self.min_call_depth:
                return "1"
            else:
                recoded = dp/float(self.min_call_depth)
                if recoded == 0.0:
                    return "0"
                else:
                    return str(recoded)

    def get_info_field_names(self):
        info_fields = OrderedDict.fromkeys(self.parser.infos.keys()).keys()
        logging.debug("(VCFRecoder) Getting following info fields:\n%s" % info_fields)
        return info_fields

    def __check_info_to_include(self):
        # Check to make sure all info columns in 'info_to_include' actually appear in VCF
        missing_columns = [x for x in self.info_to_include if x not in self.info_columns]
        if len(missing_columns) != 0:
            logging.error("INFO file error! The following INFO columns were requested but don't actually appear in the VCF:\n%s"% missing_columns)
            raise IOError("Onoe or more INFO columns specified in the info-to-include file doesn't appear in the VCF!")

    def __parse_info_to_include_file(self):
        # Parse include file to determine what information columns to include in output

        # Return if include file not specified (include all columns)
        if self.info_column_file is None:
            return None

        # Go line by line and get names of columns to include and add them to list
        info_to_include = []
        with open(self.info_column_file, "r") as info_fh:
            for line in info_fh:
                info = line.strip()
                if info != "":
                    # Check to make sure you're not just getting an empty line
                    info_to_include.append(info)

        logging.debug("Info file specifies the following columns: %s" % info_to_include)
        if len(info_to_include) == 0:
            # Return None if file didn't actually contain any data
            return None

        # Return list of info columns to include
        return info_to_include







