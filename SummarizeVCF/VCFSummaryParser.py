import logging

from VCFSummary import VCFSummary

class VCFSummaryParser:

    COUNT   = 1
    DEPTH   = 2
    QUAL    = 3
    AAFS    = 4
    INSERT  = 5
    DELETE  = 6

    def __init__(self, vcf_summary_file):
        self.vcf_summary_file = vcf_summary_file
        self.vcf_summary = None

    def parse(self):

        # Make pass to initialize VCF summary
        self.vcf_summary = self.__init_vcf_summary()

        # Make pass to fill VCF summary
        self.__fill_vcf_summary()

        # Return VCFSummary object
        return self.vcf_summary

    def __fill_vcf_summary(self):
        # Population VCFSummary with data from summary file
        curr_section = None
        next_line_is_header = False
        with open(self.vcf_summary_file, "r") as fh:
            for line in fh:
                line = line.strip()
                prev_section = curr_section
                curr_section = self.__get_section(line, curr_section)
                is_new_section = prev_section != curr_section
                if is_new_section:
                    # Skip section delimiters
                    next_line_is_header = True
                elif next_line_is_header:
                    # Skip section headers
                    next_line_is_header = False
                elif len(line) > 0:
                    # Add data for next sample (Skip empty lines)
                    sample = line.split()[0]
                    data = line.split()[1:]
                    if curr_section == VCFSummaryParser.COUNT:
                        # Set each count to the value observed in VCFSummary file
                        for i in range(len(data)):
                            count_name = self.vcf_summary.count_names[i]
                            self.vcf_summary.summary[sample].summary[count_name] = int(data[i])
                    elif curr_section == VCFSummaryParser.DEPTH:
                        # Set sample depth distribution observed in VCFSummary file
                        self.vcf_summary.summary[sample].depths = [int(x) for x in data]
                    elif curr_section == VCFSummaryParser.QUAL:
                        # Set sample qual distribution observed in VCFSummary file
                        self.vcf_summary.summary[sample].quals = [int(x) for x in data]
                    elif curr_section == VCFSummaryParser.AAFS:
                        # Set sample aaf spectrum observed in VCFSummary file
                        self.vcf_summary.summary[sample].afs = [int(x) for x in data]
                    elif curr_section == VCFSummaryParser.INSERT:
                        # Set sample insert length distribution observed in VCFSummary file
                        self.vcf_summary.summary[sample].inserts = [int(x) for x in data]
                    elif curr_section == VCFSummaryParser.DELETE:
                        # Set sample deletion length distribution observed in VCFSummary file
                        self.vcf_summary.summary[sample].deletes = [int(x) for x in data]
                    else:
                        logging.error("Invalid VCFSummary file: %s" % self.vcf_summary_file)
                        raise IOError("Invalid VCFSummary file!")

    def __init_vcf_summary(self):
        curr_section = None
        next_line_is_header = False
        header_data = None

        # VCFSummary parameters to determine
        params_to_parse = ["sample_names", "count_names", "max_depth", "max_qual", "num_afs_bins", "max_indel_len", "max_del_len"]
        summary_params = {}
        for param in params_to_parse:
            summary_params[param] = None

        with open(self.vcf_summary_file, "r") as fh:
            for line in fh:
                line = line.strip()
                prev_section = curr_section
                curr_section = self.__get_section(line, curr_section)
                is_new_section = prev_section != curr_section
                if is_new_section:
                    next_line_is_header = True
                elif next_line_is_header:
                    header_data = line.split()[1:]
                    if curr_section == VCFSummaryParser.COUNT:
                        summary_params["count_names"] = header_data
                    elif curr_section == VCFSummaryParser.DEPTH:
                        summary_params["max_depth"] = int(header_data[-2]) + 1
                    elif curr_section == VCFSummaryParser.QUAL:
                        summary_params["max_qual"] = int(header_data[-2]) + 1
                    elif curr_section == VCFSummaryParser.AAFS:
                        summary_params["num_afs_bins"] = len(header_data)
                    elif curr_section == VCFSummaryParser.INSERT:
                        summary_params["max_indel_len"] = int(header_data[-2]) + 1
                    elif curr_section == VCFSummaryParser.DELETE:
                        summary_params["max_del_len"] = int(header_data[-2]) + 1
                    elif len(line) > 0:
                        # Throw error if a non-empty line is reach and the section is unknown
                        logging.error("Invalid VCFSummary file: %s" % self.vcf_summary_file)
                        raise IOError("Invalid VCFSummary file!")

                    # Indicate that next line is not header
                    next_line_is_header = False

                # Get sample names from count section
                elif curr_section == VCFSummaryParser.COUNT:
                    sample_name = line.split()[0]
                    if summary_params["sample_names"] is None:
                        # Initialize sample name list
                        summary_params["sample_names"] = []
                    # Add sample name to list
                    summary_params["sample_names"].append(sample_name)

        # Check to make sure all params were parsed from summary file
        errors = False
        for param_name, param_value in summary_params.iteritems():
            if param_value is None:
                errors = True
                logging.error("Unable to determine parameter '%s' from VCFSummary file!" % param_name)
        if errors:
            raise IOError("Invalid VCFSummary file! One or more sections missing from input!")

        return VCFSummary(summary_params["sample_names"],
                          required_count_names=summary_params["count_names"],
                          max_depth=summary_params["max_depth"],
                          max_qual=summary_params["max_qual"],
                          num_afs_bins=summary_params["num_afs_bins"],
                          max_indel_len=summary_params["max_indel_len"])

    def __get_section(self, line, curr_section):
        # Return current section of parser based on the previous section and the current line
        if line.startswith("#COUNTS"):
            return VCFSummaryParser.COUNT
        elif line.startswith("#DEPTH"):
            return VCFSummaryParser.DEPTH
        elif line.startswith("#QUAL"):
            return VCFSummaryParser.QUAL
        elif line.startswith("#AAFS"):
            return VCFSummaryParser.AAFS
        elif line.startswith("#INSERT"):
            return VCFSummaryParser.INSERT
        elif line.startswith("#DELETE"):
            return VCFSummaryParser.DELETE
        else:
            return curr_section
