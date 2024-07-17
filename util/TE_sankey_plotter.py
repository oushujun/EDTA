import pandas as pd
import os
import numpy as np
import subprocess
import click
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import seaborn as sns
from math import pi

# Jiangzhao Qian from TEtrimmer
# Apr 2024

current_path = __file__
current_dir = os.path.dirname(os.path.abspath(current_path))


def set_te_category_and_color(less_classification=False):
    # categories are modified based on https://github.com/oushujun/EDTA/blob/master/lib-test.pl
    if not less_classification:
        categories = {
            'LTR': {"RLG", "RLC", "RLB", "RLR", "RLE", "LTR", "RLX", "Gypsy", "Copia", "BEL"},
            'LINE': {"LINE", "RIL", "RIT", "RIX", "Penelope"},
            'SINE': {"SINE", "RST", "RSX", "Non-LTR"},
            'DNA-TIR': {"DNA", "TIR", "hAT", "hAT-Ac", "MULE", "MLE", "MuDR", "Tourist", "CACTA", "PILE", "POLE",
                        "Stowaway", "TcMar-Stowaway", "PIF", "Harbinger", "Tc1", "En-Spm", "EnSpm", "CMC-EnSpm", "PiggyBac",
                        "Mirage", "P-element", "Transib", "DTA", "DTH", "DTT", "DTM", "DTC", "DTX", "DTR", "DTE", "Merlin",
                        "DTP", "DTB", "polinton", "Academ", "Crypton"},
            'MITE': {"MITE"},
            'Helitron': {"Helitron", "DHH", "DHX", "helitron"},
            #'Unknown': {"Unknown", "unknown", "Unspecified"},
            'Simple-repeat': {"Low_complexity", "repeat_region", "Simple_repeat", "Low", "Satellite"}
        }

        # All colors are colorblind friendly
        categories_color = {
            'LTR': "#33A02C",  # Dark green
            'LINE': "#E31A1C",  # Bright red
            'SINE': "#FF7F00",  # Vivid orange
            'DNA-TIR': "#1F78B4",  # Bright blue
            'MITE': "#6A3D9A",  # Deep purple
            'Helitron':  "#FB9A99",  # Soft pink
            'Unknown':  "#A6CEE3",  # Sky blue
            'Simple-repeat': "#B15928",  # Rust
            'Not-TE': "#919697"  # light grey
        }

        categories_color_lighter = {
            'LTR': "rgba(51, 160, 44, 0.4)",  # Lighter Green
            'LINE': "rgba(227, 26, 28, 0.4)",  # Lighter Red
            'SINE': "rgba(255, 127, 0, 0.4)",  # Lighter Orange
            'DNA-TIR': "rgba(31, 120, 180, 0.4)",  # Lighter Blue
            'MITE': "rgba(106, 61, 154, 0.4)",  # Lighter Purple
            'Helitron': "rgba(251, 154, 153, 0.4)",  # Lighter Pink
            'Unknown': "rgba(166, 206, 227, 0.4)",  # Lighter Sky Blue
            'Simple-repeat': "rgba(177, 89, 40, 0.4)",  # Lighter Brown
            'Not-TE': "rgba(145, 150, 151, 0.4)"  # Lighter Grey
        }

        predefined_order = ['LTR', 'LINE', 'SINE', 'DNA-TIR', 'MITE', 'Helitron', 'Unknown', 'Simple_repeat', 'Not-TE']

    else:

        categories = {
            'LTR': {"RLG", "RLC", "RLB", "RLR", "RLE", "LTR", "RLX", "Gypsy", "Copia", "BEL"},
            'Non-LTR': {"LINE", "RIL", "RIT", "RIX", "Penelope", "SINE", "RST", "RSX", "Non-LTR", "nonLTR", "YR", "Retroposon"},
            'DNA': {"DNA", "TIR", "hAT", "hAT-Ac", "MULE", "MLE", "MuDR", "Tourist", "CACTA", "PILE", "POLE",
                    "Stowaway", "TcMar-Stowaway", "PIF", "Harbinger", "Tc1", "En-Spm", "EnSpm", "CMC-EnSpm", "PiggyBac",
                    "Mirage", "P-element", "Transib", "DTA", "DTH", "DTT", "DTM", "DTC", "DTX", "DTR", "DTE", "Merlin",
                    "DTP", "DTB", "polinton", "Academ", "Crypton", "MITE", "Helitron", "DHH", "DHX", "helitron"},
            #'Unknown': {"Unknown", "unknown", "Unspecified"},
            'Simple_repeat': {"Low_complexity", "repeat_region", "Simple_repeat", "Low", "Satellite"}
        }

        categories_color = {
            'LTR': "#33A02C",  # Dark green
            'Non-LTR': "#E31A1C",  # Bright red
            'DNA': "#1F78B4",  # Bright blue
            'Unknown':  "#A6CEE3",  # Sky blue
            'Simple-repeat': "#B15928",  # Rust
            'Not-TE': "#919697"  # light grey
        }

        categories_color_lighter = {
            'LTR': "rgba(51, 160, 44, 0.4)",  # Lighter Green
            'Non-LTR': "rgba(227, 26, 28, 0.4)",  # Lighter Red
            'DNA': "rgba(31, 120, 180, 0.4)",  # Lighter Blue
            'Unknown': "rgba(166, 206, 227, 0.4)",  # Lighter Sky Blue
            'Simple-repeat': "rgba(177, 89, 40, 0.4)",  # Lighter Brown
            'Not-TE': "rgba(145, 150, 151, 0.4)"  # Lighter Grey
        }

        predefined_order = ['LTR', 'Non-LTR', 'DNA', 'Unknown', 'Simple_repeat', 'Not-TE']

        """
        # For TETrimmer H. sapiens plot
        categories = {
            'LTR': {"RLG", "RLC", "RLB", "RLR", "RLE", "LTR", "RLX", "Gypsy", "Copia", "BEL"},
            'LINE': {"LINE", "RIL", "RIT", "RIX", "Penelope"},
            'SINE': {"SINE", "RST", "RSX", "Non-LTR"},
            'DNA': {"DNA", "TIR", "hAT", "hAT-Ac", "MULE", "MLE", "MuDR", "Tourist", "CACTA", "PILE", "POLE",
                    "Stowaway", "TcMar-Stowaway", "PIF", "Harbinger", "Tc1", "En-Spm", "EnSpm", "CMC-EnSpm",
                    "PiggyBac", "Mirage", "P-element", "Transib", "DTA", "DTH", "DTT", "DTM", "DTC", "DTX", "DTR",
                    "DTE", "Merlin", "DTP", "DTB", "polinton", "Academ", "Crypton", "MITE", "Helitron", "DHH",
                    "DHX", "helitron"},
            #'DNA': {"MITE"},
            #'DNA': {"Helitron", "DHH", "DHX", "helitron"},
            #'Unknown': {"Unknown", "unknown", "Unspecified"},
            'Simple-repeat': {"Low_complexity", "repeat_region", "Simple_repeat", "Low", "Satellite"}
        }

        # All colors are colorblind friendly
        categories_color = {
            'LTR': "#33A02C",  # Dark green
            'LINE': "#E31A1C",  # Bright red
            'SINE': "#FF7F00",  # Vivid orange
            'DNA': "#1F78B4",  # Bright blue
            'MITE': "#6A3D9A",  # Deep purple
            'Helitron': "#FB9A99",  # Soft pink
            'Unknown': "#A6CEE3",  # Sky blue
            'Simple-repeat': "#B15928",  # Rust
            'Not-TE': "#919697"  # light grey
        }

        categories_color_lighter = {
            'LTR': "rgba(51, 160, 44, 0.4)",  # Lighter Green
            'LINE': "rgba(227, 26, 28, 0.4)",  # Lighter Red
            'SINE': "rgba(255, 127, 0, 0.4)",  # Lighter Orange
            'DNA': "rgba(31, 120, 180, 0.4)",  # Lighter Blue
            'MITE': "rgba(106, 61, 154, 0.4)",  # Lighter Purple
            'Helitron': "rgba(251, 154, 153, 0.4)",  # Lighter Pink
            'Unknown': "rgba(166, 206, 227, 0.4)",  # Lighter Sky Blue
            'Simple-repeat': "rgba(177, 89, 40, 0.4)",  # Lighter Brown
            'Not-TE': "rgba(145, 150, 151, 0.4)"  # Lighter Grey
        }

        predefined_order = ['LTR', 'LINE', 'SINE', 'DNA', 'MITE', 'Helitron', 'Unknown', 'Simple_repeat', 'Not-TE']
        """
    return categories, categories_color, categories_color_lighter, predefined_order


def reverse_category(categories):
    # Generate reverse_categories to be able to index the key
    reverse_categories = {}
    for category, types in categories.items():
        for type in types:
            reverse_categories[type] = category
    return reverse_categories


def format_number(num):
    if num < 1000:
        return f"{str(num)} bp"
    elif num < 1_000_000:
        return f"{num / 1_000: .1f} Kbp"
    elif num < 1_000_000_000:
        return f"{num / 1_000_000: .1f} Mbp"
    else:
        return f"{num / 1_000_000_000: .1f} Gbp"


def delete_file(*file_paths):
    for file_path in file_paths:
        if os.path.exists(file_path):
            os.remove(file_path)


def calculate_genome_length(genome_file):
    """
    Calculate the length of each sequence in a genome file in FASTA format,
    the total genome length, and write the lengths to an output file.

    :param genome_file: str, path to genome file in FASTA format
    :return: tuple(str, int), path to the output file containing sequence names and lengths, total genome length
    """
    genome_lengths = {}
    total_genome_length = 0  # Initialize total genome length
    with open(genome_file, "r") as f:
        current_seq = None
        current_length = 0
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_seq is not None:
                    genome_lengths[current_seq] = current_length
                    total_genome_length += current_length  # Update total genome length
                current_seq = line[1:].split(" ")[0]
                current_length = 0
            else:
                current_length += len(line)
        if current_seq is not None:
            genome_lengths[current_seq] = current_length
            total_genome_length += current_length  # Update total genome length for the last sequence

    # Write lengths to output file
    output_file = genome_file + ".length"
    with open(output_file, "w") as out:
        for seq_name, length in genome_lengths.items():
            out.write(f"{seq_name}\t{length}\n")
        # Optionally, write total genome length to the file or handle it differently

    return output_file, total_genome_length


def read_rm_out(input_file, output_dir, categories, reverse_categories, genome_length, debug=False):
    # Only extract desired columns
    df = pd.read_csv(input_file, sep=r'\s+', skiprows=3, usecols=[4, 5, 6, 8, 10], header=None)

    # Rename columns
    df.columns = ['chro', 'start', 'end', 'strand', 'type']

    df['length'] = df['end'] - df['start']
    # Reorder columns by specifying the order explicitly
    df = df[['chro', 'start', 'end', 'type', 'length', 'strand']]

    # Separate TE classification by '/' and match them with TE categories
    def map_type(type_string):
        # Check if input is not a string or is a non-applicable string value
        if not isinstance(type_string, str) or type_string in ['NA', 'None', '', 'N/A']:
            return 'Unknown'
        te_types = type_string.split('/')
        for t in reversed(te_types):
            if t in all_types:
                return reverse_categories[t]
        return 'Unknown'

    # Combine all TE types from the categories
    all_types = set()
    for te_category in categories.values():
        all_types = all_types.union(te_category)

    # Change C in the strand column to -
    df['strand'] = df['strand'].apply(lambda x: "-" if x == "C" else x)

    # Apply the new map_type function to each type
    df['type'] = df['type'].apply(map_type)

    # Sort dataframe by chromosome and the positions
    df_sorted = df.sort_values(by=['chro', 'start']).reset_index(drop=True)

    # Group by 'type' and sum the 'length'
    df_grouped_sum = df_sorted.groupby('type')['length'].sum().reset_index()

    # Write df to a file. Define output file name
    output_df_path = os.path.join(output_dir, f"{os.path.basename(input_file)}.bed")
    df.to_csv(output_df_path, sep="\t", index=False, header=False)

    # Run bedtools to sort file
    # Define bedtools sorted file name
    sorted_output_df_path = f"{output_df_path}_sort.bed"
    bed_sort = f"bedtools sort -i {output_df_path} -g {genome_length} > {sorted_output_df_path}"

    try:
        subprocess.run(bed_sort, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"bedtools sort error\n{e.stderr}")

    if not debug:
        delete_file(output_df_path)

    return sorted_output_df_path, df_grouped_sum


def get_complement_region(bed_file, genome_length):
    # Define output file
    bed_complement_path = os.path.join(os.path.dirname(bed_file), f"{os.path.basename(bed_file)}_com.bed")
    bed_get_complement = ["bedtools", "complement", "-i", bed_file, "-g", genome_length]

    try:
        # Execute command
        result = subprocess.run(bed_get_complement, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE,
                                text=True)

        # Write the output to the file
        with open(bed_complement_path, 'w') as f:
            f.write(result.stdout)

    except subprocess.CalledProcessError as e:
        print("bedtools complement error:\n")
        print(e.stderr)

    return bed_complement_path


def solve_overlaps(input_file_name, output_file_name):
    with open(output_file_name, 'w') as output_file:
        # Open and read the input file
        with open(input_file_name, 'r') as input_file:
            for line in input_file:
                fields = line.strip().split("\t")
                chrom = fields[0]
                types = fields[3].split(",")
                starts = list(map(int, fields[5].split(",")))
                ends = list(map(int, fields[6].split(",")))
                n = len(types)

                # Determine which intervals are completely covered by others
                to_remove = [False] * n
                for i in range(n):
                    for j in range(n):
                        if i != j and starts[j] <= starts[i] and ends[j] >= ends[i]:
                            to_remove[i] = True
                            break

                # Filtering out the marked intervals
                filtered_types = [types[i] for i in range(n) if not to_remove[i]]
                filtered_starts = [starts[i] for i in range(n) if not to_remove[i]]
                filtered_ends = [ends[i] for i in range(n) if not to_remove[i]]

                # Generate and write the processed lines to the output file
                for i in range(len(filtered_types)):
                    if i == 0:
                        output_line = f"{chrom}\t{filtered_starts[i]}\t{filtered_ends[i]}\t{filtered_types[i]}\n"
                    else:
                        output_line = f"{chrom}\t{filtered_ends[i-1] + 1}\t{filtered_ends[i]}\t{filtered_types[i]}\n"
                    output_file.write(output_line)

    return output_file_name


def bed_merge(bed_file, genome_length, debug=False):
    # Merge lines connect or close with each other. Define output file path
    output_path_no_overlap = os.path.join(os.path.dirname(bed_file), f"{os.path.basename(bed_file)}_no_overlap.bed")
    output_path_overlap = os.path.join(os.path.dirname(bed_file), f"{os.path.basename(bed_file)}_overlap.bed")
    output_path_overlap_solved = os.path.join(os.path.dirname(bed_file), f"{os.path.basename(bed_file)}_overlap_solved.bed")
    output_path_combined = os.path.join(os.path.dirname(bed_file), f"{os.path.basename(bed_file)}_combined.bed")
    output_path_combined_sort = os.path.join(os.path.dirname(bed_file), f"{os.path.basename(bed_file)}_combined_sort.bed")

    # Overlapped line looks like
    """
    scaffold_1	150336	152457	DNA,Non-LTR,Non-LTR	DNA,Non-LTR	150336,150404,151592	150433,151605,152457
    scaffold_1	157689	159647	Non-LTR,Non-LTR,LTR	LTR,Non-LTR	157689,159335,159462	159446,159467,159647
    scaffold_1	167882	168622	Non-LTR,Non-LTR,LTR	LTR,Non-LTR	167882,167899,168276	167965,168289,168622
    """
    # Not overlapped line looks like
    """
    scaffold_1	14994	16363	Non-LTR
    scaffold_1	16598	17129	LTR
    scaffold_1	20040	20458	LTR
    scaffold_1	20459	21491	LTR
    """

    # Ues the difference between the overlapped and no overlapped line to separate bed file
    bed_merge_no_overlap = f"""bedtools merge -i {bed_file} -d 0 -c 4,4,2,3 -o collapse,distinct,collapse,collapse | awk 'BEGIN{{OFS="\\t"}} $5!~/,/ {{print $1,$2,$3,$5}}' > {output_path_no_overlap}"""
    bed_merged_overlap = f"""bedtools merge -i {bed_file} -d 0 -c 4,4,2,3 -o collapse,distinct,collapse,collapse | awk 'BEGIN{{OFS="\\t"}} $5~/,/' > {output_path_overlap}"""
    combined_command_merge = [bed_merge_no_overlap, bed_merged_overlap]

    for command in combined_command_merge:
        try:
            subprocess.run(command, shell=True, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            print(f"bedtools merge error for command {command}\n"
                  f"{e.stderr}\n")
            raise Exception

    # Separate overlapped elements
    solve_overlaps(output_path_overlap, output_path_overlap_solved)

    # Define command line to concatenate overlap_solved and no_overlap file
    combine_f = f"cat {output_path_overlap_solved} {output_path_no_overlap} > {output_path_combined}"

    # Define command line to sort combined file
    sort_combined = f"bedtools sort -i {output_path_combined} -g {genome_length} > {output_path_combined_sort}"

    # Combine commands with '&&' to ensure each command runs only if the previous one succeeded
    combined_command_combine = [combine_f, sort_combined]

    for command in combined_command_combine:
        try:
            subprocess.run(command, shell=True, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            print(f"bedtools merge error for command {command}\n"
                  f"{e.stderr}\n")
            raise Exception
    if not debug:
        delete_file(output_path_no_overlap, output_path_overlap, output_path_overlap_solved, output_path_combined)

    return output_path_combined_sort


def bed_intersect(bed1, bed2, strand=False):
    # Define the output file
    bed1_name = os.path.basename(bed1)
    bed2_name = os.path.basename(bed2)
    bed_intersect_out_path = os.path.join(os.path.dirname(bed2), f"{bed1_name}_{bed2_name}_inte.bed")

    if strand:
        bed_intersect = f"bedtools intersect -a {bed1} -b {bed2} -wao -s > {bed_intersect_out_path}"
    else:
        bed_intersect = f"bedtools intersect -a {bed1} -b {bed2} -wao > {bed_intersect_out_path}"

    try:
        subprocess.run(bed_intersect, shell=True, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"bedtools intersect error\n"
              f"{e.stderr}\n")
    return bed_intersect_out_path


def parse_intersect(input_file, query_com=False, ref_com=False):
    df = pd.read_csv(input_file, sep='\t', header=None)

    # When query and reference are both not complementary bed file
    if not query_com and not ref_com:
        # Filter out rows where columns 7 and 8 are equal to -1
        filtered_df = df[(df[5] != -1) & (df[6] != -1)]
        # Group by column 4 and column 10, then sum the last column
        grouped_sum = filtered_df.groupby([3, 7])[len(df.columns) - 1].sum().reset_index()

    # When only reference is complementary bed file
    elif not query_com and ref_com:
        # Filter out rows where columns 7 and 8 are equal to -1
        filtered_df = df[(df[5] != -1) & (df[6] != -1)]
        grouped_sum = filtered_df.groupby(3)[len(df.columns) - 1].sum().reset_index()
        grouped_sum.insert(1, 'reference', 'Not-TE')  # Inserting the 'reference' column

    # When only query is complementary bed file
    elif query_com and not ref_com:
        # Filter out rows where columns 4 and 5 are equal to -1
        filtered_df = df[(df[4] != -1) & (df[5] != -1)]
        grouped_sum = filtered_df.groupby(6)[len(df.columns) - 1].sum().reset_index()
        grouped_sum.insert(0, 'query', 'Not-TE')  # Inserting the 'query' column

    # When query and reference are both complementary bed file
    else:
        # Filter out rows where columns 4 and 5 are equal to -1
        filtered_df = df[(df[4] != -1) & (df[5] != -1)]
        grouped_sum_n = filtered_df[len(df.columns) - 1].sum()

        # Create dataframe. Because grouped_sum_n is only a number
        grouped_sum = pd.DataFrame({
            'query': ['Not-TE'],
            'reference': ['Not-TE'],
            'length': [grouped_sum_n]
        })

    grouped_sum.columns = ['query', 'reference', 'length']

    return grouped_sum


def bed_merge_and_cross_intersect(query_bed, reference_bed, genome_length, debug=False):
    #####################################################################################################
    # Code block: Merge and get complementary bed files
    #####################################################################################################

    # Merge bed file for query
    merged_query_bed = bed_merge(query_bed, genome_length, debug=debug)
    # Create complement bed file for query
    query_bed_complement = get_complement_region(merged_query_bed, genome_length)
    # Merge bed file for reference
    merged_reference_bed = bed_merge(reference_bed, genome_length, debug=debug)
    # Create complement bed file for reference
    reference_bed_complement = get_complement_region(merged_reference_bed, genome_length)

    #####################################################################################################
    # Code block: bed file interactions
    #####################################################################################################

    # Do bed intersection for query with reference
    query_ref_inter = bed_intersect(merged_query_bed, merged_reference_bed, strand=False)
    # Do bed intersection for query with reference complement
    query_ref_com_inter = bed_intersect(merged_query_bed, reference_bed_complement, strand=False)
    # Do bed intersection for query complement with reference
    query_com_ref_inter = bed_intersect(query_bed_complement, merged_reference_bed, strand=False)
    # Do bed intersection for query complement with reference complement
    query_com_ref_com_inter = bed_intersect(query_bed_complement, reference_bed_complement, strand=False)

    #####################################################################################################
    # Code block: parse intersections
    #####################################################################################################

    # Calculate sum
    query_ref_inter_df = parse_intersect(query_ref_inter, query_com=False, ref_com=False)
    query_ref_com_inter_df = parse_intersect(query_ref_com_inter, query_com=False, ref_com=True)
    query_com_ref_inter_df = parse_intersect(query_com_ref_inter, query_com=True, ref_com=False)
    query_com_ref_com_inter_df = parse_intersect(query_com_ref_com_inter, query_com=True, ref_com=True)

    # Combine DataFrames and sort by query
    combined_df = pd.concat([query_ref_inter_df, query_ref_com_inter_df, query_com_ref_inter_df,
                             query_com_ref_com_inter_df], ignore_index=True)

    combined_df['query'] = combined_df['query'].astype(str)
    combined_df = combined_df.sort_values(by='query', ignore_index=True)

    # Delete file
    if not debug:
        delete_file(merged_query_bed, query_bed_complement, merged_reference_bed, reference_bed_complement,
                    query_ref_inter, query_ref_com_inter, query_com_ref_inter, query_com_ref_com_inter)

    return combined_df


def sankey_plot(combined_df, categories_color, categories_color_lighter, output_path, show_plot=False):
    length_query = combined_df.groupby('query')['length'].sum()

    length_reference = combined_df.groupby('reference')['length'].sum()

    # Get unique labels for query and reference
    query_labels = combined_df['query'].unique()

    query_labels_with_length = []
    for i in query_labels:
        type_l = format_number(length_query[i])
        type_with_l = f"{i} {type_l}"
        query_labels_with_length.append(type_with_l)

    query_labels_dict = {label: idx for idx, label in enumerate(query_labels)}
    reference_labels = combined_df['reference'].unique()

    reference_labels_with_length = []
    for i in reference_labels:
        type_l = format_number(length_reference[i])
        type_with_l = f"{i} {type_l}"
        reference_labels_with_length.append(type_with_l)

    reference_labels_dict = {label: idx for idx, label in enumerate(reference_labels, start=len(query_labels))}

    all_labels = np.append(query_labels, reference_labels)
    all_labels_with_length = np.append(query_labels_with_length, reference_labels_with_length)

    # Generate source, target, and value lists for the Sankey diagram
    source = combined_df['query'].map(query_labels_dict)
    target = combined_df['reference'].map(reference_labels_dict)
    values = combined_df['length']

    # Map colors for each node
    node_colors = [categories_color.get(label, "#000000") for label in all_labels]  # Default to black if not found

    # Map target indices to their corresponding node colors for links
    node_colors_lighter = [categories_color_lighter.get(label, "#000000") for label in all_labels]
    link_colors = [node_colors_lighter[idx] for idx in target]

    # Create the Sankey diagram
    fig = go.Figure(data=[go.Sankey(
        arrangement='snap',
        node=dict(
            pad=10,
            thickness=20,
            line=dict(color="black", width=2),
            label=all_labels_with_length,
            color=node_colors
        ),
        link=dict(
            source=source,
            target=target,
            value=values,
            color=link_colors
        )
    )])

    fig.update_layout(title_text= "Left side represents query TE annotation and right side is standard TE annotation.",
                      font_size=13)
    if show_plot:
        fig.show()
    # Save the plot as an HTML file
    fig.write_html(output_path)


def sankey_plot_for_two_df(combined_df1, combined_df2, categories_color, categories_color_lighter):
    # Deprecated don't use this
    # Get unique labels for query and reference from combined_df1
    query_labels_df1 = combined_df1['query'].unique()
    reference_labels_df1 = combined_df1['reference'].unique()

    # Get unique labels for query and reference from combined_df2
    query_labels_df2 = combined_df2['query'].unique()
    reference_labels_df2 = combined_df2['reference'].unique()

    query_labels_df1_dict = {label: idx for idx, label in enumerate(query_labels_df1)}

    # Combine reference_labels_df1 and query_labels_df2
    reference_labels_df1_query_labels_df2 = np.append(reference_labels_df1, query_labels_df2)
    reference_labels_df1_query_labels_df2 = np.unique(reference_labels_df1_query_labels_df2)


    reference_labels_df1_query_labels_df2_dict = {label: idx for idx, label in
                                                  enumerate(reference_labels_df1_query_labels_df2,
                                                            start=len(query_labels_df1))}

    reference_labels_df2_dict = {label: idx for idx, label in
                                 enumerate(reference_labels_df2,
                                           start=(len(query_labels_df1) + len(reference_labels_df1_query_labels_df2)))}

    all_labels = np.concatenate([query_labels_df1, reference_labels_df1_query_labels_df2, reference_labels_df2])

    # Generate source, target, and value lists for the Sankey diagram
    source = np.append(combined_df1['query'].map(query_labels_df1_dict), combined_df2['query'].map(reference_labels_df1_query_labels_df2_dict))
    target = np.append(combined_df1['reference'].map(reference_labels_df1_query_labels_df2_dict), combined_df2['reference'].map(reference_labels_df2_dict))
    values = np.append(combined_df1['length'], combined_df2['length'])

    # Map colors for each node
    node_colors = [categories_color.get(label, "#000000") for label in all_labels]  # Default to black if not found

    # Map target indices to their corresponding node colors for links
    node_colors_lighter = [categories_color_lighter.get(label, "#000000") for label in all_labels]
    link_colors = [node_colors_lighter[idx] for idx in target]

    # Create the Sankey diagram
    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=10,
            thickness=20,
            line=dict(color="black", width=2),
            label=all_labels,
            color=node_colors
        ),
        link=dict(
            source=source,
            target=target,
            value=values,
            color=link_colors
        )
    )])

    fig.update_layout(title_text="Left side represents query TE annotation and right side is standard TE annotation.",
                      font_size=11)

    fig.show()
    return fig


def matrix_plot(combined_df, matrix_path):
    pivot_df = combined_df.pivot(index='query', columns='reference', values='length').fillna(0)

    # Normalize by each row (query)
    row_totals = pivot_df.sum(axis=0)
    pivot_df_normalized_by_query = pivot_df.div(row_totals, axis=1)

    # Plot heatmap
    plt.figure(figsize=(12, 10))
    sns.heatmap(pivot_df_normalized_by_query, annot=True, fmt='.2%', cmap='viridis')
    plt.xlabel('Reference')
    plt.ylabel('Query')
    plt.title('Proportion of each query TE types within each standard TE types')
    # Save the figure
    plt.savefig(matrix_path, dpi=300, bbox_inches='tight')
    plt.close()


def calculate_confusion_metrics(df):
    metrics_single = {}
    metrics = {}
    query_types = df['query'].unique()

    for query_type in query_types:
        tp = df[(df['query'] == query_type) & (df['query'] == df['reference'])]['length'].sum()
        fp = df[(df['query'] == query_type) & (df['query'] != df['reference'])]['length'].sum()
        fn = df[(df['reference'] == query_type) & (df['query'] != df['reference'])]['length'].sum()
        tn = df[(df['query'] != query_type) & (df['reference'] != query_type)]['length'].sum()

        # Store parameter to dictionary
        metrics_single[query_type] = {
            "True_positive": tp,
            "False_positive": fp,
            "False_negative": fn,
            "True_negative": tn
        }

        sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
        specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
        accuracy = (tp + tn) / (tp + fp + fn + tn) if (tp + fp + fn + tn) > 0 else 0
        precision = tp / (tp + fp) if (tp + fp) > 0 else 0
        f1_score = 2 * (precision * sensitivity) / (precision + sensitivity) if (precision + sensitivity) > 0 else 0
        fdr = 1 - precision

        # Add calculated metrics for each query type to the dictionary
        metrics[query_type] = {
            'Sensitivity': sensitivity,
            'Specificity': specificity,
            'Accuracy': accuracy,
            'Precision': precision,
            'F1_Score': f1_score,
            'FDR': fdr
        }

    # Convert the dictionary to a DataFrame after collecting all metrics
    metrics_single_df = pd.DataFrame(metrics_single).T
    metrics_df = pd.DataFrame(metrics).T  # Transpose to have query types as the index and metrics as columns

    return metrics_df, metrics_single_df


def calculate_all_te_confusion_metrics(df):
    # This function will combine all TE types as TE and calculate the confusion metrics
    # Change all not Not-TE label to TE
    df['query_reclassified'] = np.where(df['query'] == 'Not-TE', 'Not-TE', 'All-TE')
    df['reference_reclassified'] = np.where(df['reference'] == 'Not-TE', 'Not-TE', 'All-TE')

    # Initialize metrics dictionary
    metrics_parameter = {}
    metrics = {}

    # Calculate metrics for TE and Not-TE
    for query_type in ['All-TE']:
        tp = df[(df['query_reclassified'] == query_type) & (df['query_reclassified'] == df['reference_reclassified'])]['length'].sum()
        fp = df[(df['query_reclassified'] == query_type) & (df['query_reclassified'] != df['reference_reclassified'])]['length'].sum()
        fn = df[(df['reference_reclassified'] == query_type) & (df['query_reclassified'] != df['reference_reclassified'])]['length'].sum()
        tn = df[(df['query_reclassified'] != query_type) & (df['reference_reclassified'] != query_type)]['length'].sum()

        sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
        specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
        accuracy = (tp + tn) / (tp + fp + fn + tn) if (tp + fp + fn + tn) > 0 else 0
        precision = tp / (tp + fp) if (tp + fp) > 0 else 0
        f1_score = 2 * (precision * sensitivity) / (precision + sensitivity) if (precision + sensitivity) > 0 else 0
        fdr = 1 - precision

        # Store parameter to dictionary
        metrics_parameter[query_type] = {
            "True_positive": tp,
            "False_positive": fp,
            "False_negative": fn,
            "True_negative": tn
        }

        # Add calculated metrics for each query type to the dictionary
        metrics[query_type] = {
            'Sensitivity': sensitivity,
            'Specificity': specificity,
            'Accuracy': accuracy,
            'Precision': precision,
            'F1_Score': f1_score,
            'FDR': fdr
        }

    # Convert the metrics dictionary to a DataFrame
    metrics_parameter_df = pd.DataFrame(metrics_parameter).T
    metrics_df = pd.DataFrame(metrics).T

    return metrics_df, metrics_parameter_df


def spider_plot(metrics_df, output_path, categories_color, predefined_order):
    # Set 'type' as the index
    metrics_df.set_index('type', inplace=True)

    # Reorder the dataframe according to predefined_order
    metrics_df = metrics_df.reindex(predefined_order)

    # Reset index to get 'type' back as a column
    metrics_df.reset_index(inplace=True)

    # Number of variables
    categories = list(metrics_df)[1:]
    N = len(categories)

    # What will be the angle of each axis in the plot?
    angles = [n / float(N) * 2 * pi for n in range(N)]
    angles += angles[:1]

    # Define subplot grid layout
    n_types = metrics_df.shape[0]
    n_cols = 4  # Adjust as needed
    n_rows = int(np.ceil(n_types / n_cols))

    fig, axs = plt.subplots(n_rows, n_cols, figsize=(n_cols * 3, n_rows * 3), subplot_kw=dict(polar=True))

    # Flatten axes array if it is multidimensional
    axs = axs.flatten()

    for i, row in metrics_df.iterrows():
        ax = axs[i]
        values = row.drop('type').values.flatten().tolist()
        values += values[:1]
        color = categories_color.get(row['type'], 'black')  # Fallback color

        ax.plot(angles, values, linewidth=2, linestyle='solid', label=row['type'], color=color)
        ax.fill(angles, values, color=color, alpha=0.1)
        ax.set_title(row['type'], size=13, color=color, position=(0.5, 1.1))
        ax.set_ylim(0, 1)

        # Set xticks and labels for each subplot
        ax.set_xticks(angles[:-1], categories, color='grey', size=7)
        ax.tick_params(axis='x', pad=12)
        ax.set_yticks([0.2, 0.4, 0.6, 0.8])
        ax.set_yticklabels(["0.2", "0.4", "0.6", "0.8"], color="grey", size=7)

    # Hide unused subplots
    for i in range(n_types, n_rows * n_cols):
        axs[i].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    #plt.show()


@click.command(context_settings=dict(max_content_width=120),
               help="""\b
                ##########################################################################################
                Transposable element Sankey plotter and confusion metrics calculator.
                  
                Version: v1.0.1 (21/Mar/2024) 

                Github: https://github.com/qjiangzhao/TETrimmer

                Developers:                                                                                                       
                Jiangzhao Qian;      RWTH Aachen University;             Email: jqian@bio1.rwth-aachen.de                                                                                                 
                ##########################################################################################              

""")
@click.option("--query", "-q", required="True", type=str,
              help="Query RepeatMasker out file. The out file you want to use to compare with the standard.")
@click.option("--reference", "-r", required="True", type=str,
              help="RepeatMasker out file generated by standard TE library (manually curated TE library).")
@click.option("--output_dir", "-o", default=None, type=str,
              help="Output directory. By default, it is query file path.")
@click.option("--genome_dir", "-g", default=None, type=str,
              help="Genome file path. You have to supply genome file or genome length file. If only genome "
                   "file path is supplied, the genome length can be calculated automatically.")
@click.option("--genome_length", "-gl", default=None, type=str,
              help="Genome length file path. File should not contain header, the first column is 'Chromosome' "
                   "and second column is 'length'. They are separated by tab.")
@click.option('--less_classification', "-lc", default=False, is_flag=True,
              help='Use this option to have less TE classification in the plot. For example, merge "LINE" and "SINE"'
                   'to "Non-LTR.')
@click.option('--debug', default=False, is_flag=True,
              help='Use to enable debug mode. This will keep all raw files.')
def calculate_overlap(query, reference, output_dir, genome_dir, genome_length, less_classification=False,
                      debug=False, query2=None):
    # com: complementary
    # me: merge
    # inte: intersect

    click.echo("TE sankey plotter is running......")

    if output_dir is None:
        output_dir = os.path.dirname(query)
    os.makedirs(output_dir, exist_ok=True)

    # Get genome length file and total genome length
    if genome_dir is not None:
        genome_length_f, total_genome_length = calculate_genome_length(genome_dir)
    else:
        genome_length_f = genome_length

    # Create output directory if output directory is not exist
    os.makedirs(output_dir, exist_ok=True)

    categories, categories_color, categories_color_lighter, predefined_order = set_te_category_and_color(less_classification)

    # Set reverse categories
    reverse_categories = reverse_category(categories)

    # Define sankey plot output path
    sankey_path = os.path.join(output_dir, f"{os.path.basename(query)}_sankey.html")

    # Define matrix path
    matrix_path = os.path.join(output_dir, f"{os.path.basename(query)}_heatmap_proportion_plot.pdf")

    # Define combined dataframe path
    combined_df_path = os.path.join(output_dir, f"{os.path.basename(query)}_TE_length.txt")

    # Define confusion metrics path
    confusion_metrics_path = os.path.join(output_dir, f"{os.path.basename(query)}_confusion_metrics.txt")

    # Define confusion metrics parameter path
    confusion_metrics_path_par = os.path.join(output_dir, f"{os.path.basename(query)}_confusion_metrics_TP_TN_FP_FN.txt")

    # Define spider plot output path
    spider_out_path = os.path.join(output_dir, f"{os.path.basename(query)}_spider_plot.pdf")

    # Unify TE types and generate sorted bed file
    query_bed, query_sum = read_rm_out(query, output_dir, categories, reverse_categories, genome_length_f, debug=debug)
    reference_bed, reference_sum = read_rm_out(reference, output_dir, categories, reverse_categories, genome_length_f,
                                               debug=debug)

    # Deprecated don't use this
    if query2 is not None:
        query2_bed, query2_sum = read_rm_out(query2, output_dir, categories, reverse_categories, genome_length_f)
        combined_df_query1_query2 = bed_merge_and_cross_intersect(query_bed, query2_bed, genome_length_f)
        combined_df_query2_ref = bed_merge_and_cross_intersect(query2_bed, reference_bed, genome_length_f)
        fig = sankey_plot_for_two_df(combined_df_query1_query2, combined_df_query2_ref, categories_color,
                                     categories_color_lighter)
        # Define output figure name
        fig_output_path = os.path.join(output_dir, f"fig_test.pdf")
        fig.write_image(fig_output_path)

    else:
        combined_df = bed_merge_and_cross_intersect(query_bed, reference_bed, genome_length_f, debug=debug)

        # Write combined_df to file
        combined_df.to_csv(combined_df_path, sep="\t", index=False)

        separated_te_confusion_metrics, separated_te_confusion_metrics_par = calculate_confusion_metrics(combined_df)
        confusion_metrics_all_te, confusion_metrics_all_te_par = calculate_all_te_confusion_metrics(combined_df)

        # Combine separated_te_confusion_metrics and confusion_metrics_all_te
        combined_confusion_metrics = pd.concat([separated_te_confusion_metrics, confusion_metrics_all_te])
        combined_confusion_metrics_par = pd.concat([separated_te_confusion_metrics_par, confusion_metrics_all_te_par])

        # Convert index to column
        # Reset the index of the DataFrame, moving the index to a column
        combined_confusion_metrics_reset = combined_confusion_metrics.reset_index()
        combined_confusion_metrics_reset = combined_confusion_metrics_reset.rename(columns={'index': 'type'})

        combined_confusion_metrics_par_reset = combined_confusion_metrics_par.reset_index()
        combined_confusion_metrics_par_reset = combined_confusion_metrics_par_reset.rename(columns={'index': 'type'})

        # Write combined_confusion_metrics to a file
        combined_confusion_metrics_reset.to_csv(confusion_metrics_path, sep="\t", index=False)
        combined_confusion_metrics_par_reset.to_csv(confusion_metrics_path_par, sep="\t", index=False)

        # Do sankey plot
        sankey_plot(combined_df, categories_color, categories_color_lighter, sankey_path)

    # Do pivot plot
    matrix_plot(combined_df, matrix_path)

    # Do spider plot
    spider_predefined_order = predefined_order.append('All_TE')
    spider_plot(combined_confusion_metrics_reset, spider_out_path, categories_color, spider_predefined_order)

    if not debug:
        delete_file(query_bed, reference_bed)

    click.echo(f"TE sankey plotter is finished.")


if __name__ == '__main__':
    calculate_overlap()
