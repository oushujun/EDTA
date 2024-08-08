#TODO only include one target site duplication for higher level TEs
#TODO dealwith the case that one TE contain two inserted TEs
#TODO define TE direction for those repeat region with ? at strand column, I found this is
#       not necessary because cd-hit-est will consider direction inforamtion. Do domain prediction
#       for all consensus sequence at the end will be sufficient
import re

class EDTAGFFProcessor:
    """
    Class will clean EDTA intact gff file. Separate nested TEs to single elements.
    """
    def __init__(self, input_file_path, output_file_path):
        # Initialize instance variables for input and output file paths
        self.input_file_path = input_file_path
        self.output_file_path = output_file_path
        self.gff_data = []
        self.remove_list = []
        self.result = [] # a list contain filtered gff file (only nested)
        self.result_filter = []
        self.read_gff_data()
        self.process_gff_data()
        self.remove_solo_ltr()
        self.update_direction()

    def read_gff_data(self):
        # Read the GFF3 file and return its contents as a list of lines
        # Each lines will be separated and stored into gff_file as a list
        with open(self.input_file_path, 'r') as gff_file:
            self.gff_data = gff_file.read().split('\n')
        print("finish read_EDTA_intact_gff_data")

    def write_gff_data(self):
        # Write the processed GFF3 data to the output file
        with open(self.output_file_path, 'w') as gff_file:
            gff_file.write('\n'.join(self.result_filter))
        print("finish write_EDTA_intact_gff_data")

    def process_gff_data(self):
        # Process the GFF3 data by searching for nested LTR_retrotransposon sequences
        # and splitting the outer LTR_retrotransposon into two parts accordingly

        for idx, line in enumerate(self.gff_data):
            fields = line.split() # split each gff line by white space
            if not fields: # Check if fields is empty
                continue
            chrom = fields[0] # Store gff features to variables
            seq_type = fields[2]
            start = int(fields[3])
            end = int(fields[4])
            LTR_ID = fields[8].split(";")[0]  # split gff column 9 by ";" and store features to variables
            remove_ID = fields[8].split(";")[1]
            # we only concern nested LTR-RT
            if "LTR_retrotransposon" in seq_type:
                nested_te = None
                for inner_line in self.gff_data: # start a new loop to seek nested LTR-RT
                    inner_fields = inner_line.split()
                    if not inner_fields:
                        continue
                    inner_chrom = inner_fields[0]
                    inner_seq_type = inner_fields[2]
                    inner_start = int(inner_fields[3])
                    inner_end = int(inner_fields[4])
                    # according to start and end position to see nested situation
                    if "LTR_retrotransposon" in inner_seq_type and inner_start > start and inner_end < end and inner_chrom == chrom:
                        nested_te = (inner_start, inner_end)
                        break
                    # filter out TEs that share same LTR and eliminate longer one
                    elif "LTR_retrotransposon" in inner_seq_type and inner_start == start and inner_end < end and inner_chrom == chrom:
                        self.remove_list.append(remove_ID.replace("Parent=repeat_", ""))
                        break
                    elif "LTR_retrotransposon" in inner_seq_type and inner_start > start and inner_end == end and inner_chrom == chrom:
                        self.remove_list.append(remove_ID.replace("Parent=repeat_", ""))
                        break

                if nested_te: # initially nested_te is None, once it contains a vale it will trim nested TE
                    self.result.append(
                        re.sub(f"{end}", f"{nested_te[0] - 1}", line).replace(f"{LTR_ID}", f"{LTR_ID}__1"))
                    self.result.append(
                        re.sub(f"{start}", f"{nested_te[1] + 1}", line).replace(f"{LTR_ID}", f"{LTR_ID}__2"))
                else:
                    self.result.append(line)
            else: # for lines don't contain LTR_retrotransposon, simply append it the result list.
                self.result.append(line)
        print("finish process")

    # remove according remove_list
    def remove_solo_ltr(self):
        self.result_filter = [item for item in self.result if not any(pattern in item for pattern in self.remove_list)]

    # convert ? as + in the gff file.
    # it isn't neceaasry to determine the right direction now. cd-hit-est will use both direction for clustering.
    # after MSA, use BLAST module determine all protein domains and find the right direction
    def update_direction(self):
        # Update the direction column for lines marked with "?"
        for idx, line in enumerate(self.result_filter):
            fields = line.split()
            if fields[6] == "?":
                fields[6] = "+"
                line = "\t".join(fields)
                self.result_filter[idx] = line


