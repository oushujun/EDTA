#!/usr/bin/env python3
# Tianyu Lu (tlu83@wisc.edu)
# 2023-09-22

import argparse
import os
import shutil
import tempfile

import sys
sys.path.insert(0, f"{os.path.dirname(__file__)}/bin")
from bin.main import TIRLearner


if __name__ == "__main__":

    # ================================================ argument parsing ================================================
    parser = argparse.ArgumentParser(prog="TIR-Learner")
    parser.add_argument("-f", "--genome_file", help="Genome file in fasta format", required=True)
    parser.add_argument("-n", "--genome_name", help="Genome name (Optional)", default="TIR-Learner")
    parser.add_argument("-s", "--species", help="One of the following: \"maize\", \"rice\" or \"others\"",
                        required=True)
    parser.add_argument("-g", "--grfp", help="Path to GRF program (Optional)",
                        default=os.path.dirname(shutil.which("grf-main")))
    parser.add_argument('-m', '--mode', help=("Execution mode of GRF, one of the following: \"boost\", \"strict\""
                                              " or \"Native\" (Optional, default is \"boost\")"), default="boost")
    parser.add_argument("-l", "--length", help="Max length of TIR (Optional)", default=5000)
    parser.add_argument("-t", "--processor", help="Number of processors allowed (Optional)", default=os.cpu_count())
    parser.add_argument("-o", "--output_dir", help="Output directory (Optional)", default=None)
    parser.add_argument('-d', '--debug', help="Output each module's result in csv file (Optional)", action="store_true")
    parser.add_argument('-v', '--verbose', help="Verbose mode, will show interactive progress bar (Optional)",
                        action="store_true")
    parsed_args = parser.parse_args()

    genome_file = parsed_args.genome_file
    genome_name = parsed_args.genome_name
    species = parsed_args.species
    GRF_path = parsed_args.grfp
    GRF_mode = parsed_args.mode

    TIR_length = int(parsed_args.length)
    cpu_cores = int(parsed_args.processor)
    output_dir = parsed_args.output_dir
    if output_dir is None:
        output_dir = os.path.dirname(genome_file)

    flag_debug = parsed_args.debug
    flag_verbose = parsed_args.verbose

    # Transforming the possible relative path into absolute path
    genome_file = os.path.abspath(str(genome_file))
    output_dir = os.path.abspath(str(output_dir))
    # ==================================================================================================================

    temp_dir = tempfile.mkdtemp()
    genome_file_soft_link = os.path.join(temp_dir, "genome_file_soft_link.fa.lnk")
    os.symlink(genome_file, genome_file_soft_link)
    os.chdir(temp_dir)

    TIRLearner_instance = TIRLearner(genome_file, genome_name, output_dir, species, TIR_length,
                                     GRF_path, cpu_cores, GRF_mode, flag_verbose, flag_debug)

    shutil.rmtree(temp_dir)
