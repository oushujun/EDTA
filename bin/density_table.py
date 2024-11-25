import argparse
import subprocess
import os
import sys
import shutil

def create_output_directory(genome_file):
    """
    Checks if the output directory based on the genome file name already exists.
    If it does, exits the script with an error message.
    """
    directory_name = genome_file.rsplit('.', 1)[0] + "_density_table"
    if os.path.exists(directory_name):
        sys.exit("ERROR: Intermediate directory '{}' already exists. Please remove it or run the program from a different path.".format(directory_name))
    else:
        os.makedirs(directory_name)
    return directory_name

def reformat_gff(input_gff, output_file, output_directory):
    """
    Reformats the GFF file to include only a subset of columns.
    """
    with open(input_gff, 'r') as infile, open(os.path.join(output_directory, output_file), 'w') as outfile:
        for line in infile:
            if not line.startswith('#'):
                parts = line.strip().split('\t')
                outfile.write(f"{parts[0]}\t{parts[3]}\t{parts[4]}\t{parts[2]}\n")

def create_bed_files(input_file, output_directory):
    """
    Creates .bed files for each type found in the GFF file.
    """
    with open(os.path.join(output_directory, input_file), 'r') as infile:
        for line in infile:
            parts = line.strip().split('\t')
            with open(os.path.join(output_directory, parts[3] + ".bed"), 'a') as outfile:
                outfile.write(line)

def create_genome_fai(genome_file, output_directory):
    """
    Creates a genome index file using samtools.
    """
    subprocess.run(["samtools", "faidx", genome_file], stderr=subprocess.DEVNULL)
    # Move the .fai file to the output directory
    fai_file = os.path.abspath(genome_file) + ".fai"
    if os.path.exists(fai_file):
        shutil.move(fai_file, output_directory)

def make_windows(fai_file, output_file, output_directory):
    """
    Generates window sizes across the genome using bedtools.
    """
    modified_fai = "genome.fa.fai2"
    with open(fai_file, 'r') as infile, open(os.path.join(output_directory, modified_fai), 'w') as outfile:
        for line in infile:
            parts = line.split('\t')
            outfile.write(f"{parts[0]}\t{parts[1]}\n")

    with open(os.path.join(output_directory, output_file), 'w') as outfile:
        subprocess.run(["bedtools", "makewindows", "-g", modified_fai, "-w", "1000000", "-s", "500000"], cwd=output_directory, stdout=outfile, stderr=subprocess.DEVNULL)

def calculate_density(window_file, output_directory):
    """
    Calculates the density for each TE type using bedtools coverage.
    """
    for bed_file in os.listdir(output_directory):
        if bed_file.endswith('.bed') and bed_file != window_file:
            type_name = bed_file.replace('.bed', '')
            density_file = type_name + '.density.bed'
            with open(os.path.join(output_directory, density_file), 'w') as outfile:
                subprocess.run(["bedtools", "coverage", "-a", window_file, "-b", bed_file], cwd=output_directory, stdout=outfile, stderr=subprocess.DEVNULL)

def merge_and_format(output_directory):
    """
    Merges and formats the output.
    """
    for density_file in os.listdir(output_directory):
        if density_file.endswith('.density.bed'):
            type_name = density_file.replace('.density.bed', '')
            with open(os.path.join(output_directory, density_file), 'r') as infile:
                for line in infile:
                    parts = line.strip().split('\t')
                    output_line = f"{parts[0]}\t{parts[1]}\t{parts[6]}\t{type_name}"
                    print(output_line)

def cleanup_files(output_directory, file_list):
    """
    Deletes specified files in the output directory.
    """
    for file_name in file_list:
        file_path = os.path.join(output_directory, file_name)
        try:
            if os.path.isfile(file_path):
                os.remove(file_path)
        except Exception as e:
            print(f"Error encountered while deleting file {file_name}: {e}", file=sys.stderr)

def delete_output_directory(directory):
    """
    Deletes the specified output directory and all its contents.
    """
    shutil.rmtree(directory)

def main():
    parser_description = (
        "Extracts density information for each TE superfamily from the EDTA GFF files.\n"
        "Requires python3, bedtools, and samtools. All can be installed with conda.\n"
        "Be sure to run this in a directory where there are not already bed files present."
    )
    parser = argparse.ArgumentParser(description=parser_description)
    parser.add_argument('-genome', type=str, required=True, help='Path to the genome file (genome.fa)')
    parser.add_argument('-gff', type=str, required=True, help='Path to the EDTA GFF file (genome.mod.EDTA.TEanno.gff3)')

    args = parser.parse_args()

    output_directory = create_output_directory(args.genome)

    # List of intermediate files to delete
    intermediate_files = [
        "reformat.gff",
        os.path.basename(args.genome) + ".fai",
        "genome.fa.fai2",
        "windows.bed.sizes"
    ]

    try:
        reformatted_gff = "reformat.gff"
        reformat_gff(args.gff, reformatted_gff, output_directory)
        create_bed_files(reformatted_gff, output_directory)
        create_genome_fai(args.genome, output_directory)
        window_file = "windows.bed.sizes"
        fai_file_path = os.path.join(output_directory, os.path.basename(args.genome) + ".fai")
        make_windows(fai_file_path, window_file, output_directory)
        calculate_density(window_file, output_directory)
        merge_and_format(output_directory)
    except Exception as e:
        print(f"Error encountered: {e}", file=sys.stderr)
        parser.print_help()
        sys.exit(1)
    finally:
        # Perform cleanup at the end
        cleanup_files(output_directory, intermediate_files)
        # Delete the output directory
        delete_output_directory(output_directory)

if __name__ == "__main__":
    main()
