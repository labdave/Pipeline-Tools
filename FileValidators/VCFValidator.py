import logging

class VCFValidator:
    # Class for validating VCF file structure

    # Canoncial columns in the order they're supposed to appear in a VCF
    HEADER_FIELDS = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]

    @staticmethod
    def is_valid(path):

        num_samples = 0

        # Check to make sure VCF file is valid format
        with open(path, "r") as vcf_fh:

            # Boolean state variable determining which checks have been performed
            state = "FormatDeclaration"

            for line in vcf_fh:

                if state == "FormatDeclaration":

                    # Check to make sure the first line says ##fileformat=VCFv4.X
                    if not line.startswith("##fileformat=VCF"):
                        logging.error("Invalid VCF file format! First line of file doesn't look like '##fileformat=VCF'")
                        return False
                    state = "MetadataHeader"

                elif state == "MetadataHeader":

                    # Check format of metadata header section
                    if line.startswith("#CHROM"):
                        # Check to see that required fields are present in header line
                        line = line.lstrip("#")
                        line = line.split()

                        # Get the number of samples in the VCF
                        num_samples = len(line) - len(VCFValidator.HEADER_FIELDS)

                        expected_header_string  = "".join(VCFValidator.HEADER_FIELDS)
                        actual_header_string    = "".join(line[0:len(VCFValidator.HEADER_FIELDS)])

                        # Check to make sure the header fields are present in the correct order
                        if actual_header_string != expected_header_string:
                            logging.error("Invalid VCF file format! Invalid static column lables.\n"
                                          "Expected: %s\n"
                                          "Received: %s" % (expected_header_string, actual_header_string))
                            return False

                        # Update current state
                        state = "FirstRecord"

                    elif not line.startswith("##"):
                        # Return false if header line doesn't start with '##'
                        logging.error("Invalid VCF file format! Metadata header lines must begin with '##'!")
                        return False

                elif state == "FirstRecord":

                    # Validate first record
                    line = line.split()
                    if len(line) != num_samples + len(VCFValidator.HEADER_FIELDS):
                        logging.error("Invalid VCF file format! First record did not contain the correct number of columns!")
                        return False

                    # Return True as VCF has been validated to the point that I give a shit about
                    logging.debug("VCF file is valid: %s" % path)
                    return True