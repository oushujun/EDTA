#!/usr/bin/env python3
# Tianyu Lu (tlu83@wisc.edu)
# 2023-09-22

import argparse
import os
import shutil

import sys
sys.path.insert(0, f"{os.path.dirname(__file__)}/bin")

# Use if True to suppress the PEP8: E402 warning
if True:  # noqa: E402
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
    parser.add_argument('-m', '--mode', help=("Execution mode of GRF, one of the following: \"boost\", \"mix\""
                                              " or \"native\" (Optional)"), default="smart")
    parser.add_argument("-l", "--length", help="Max length of TIR (Optional)", default=5000)
    parser.add_argument("-t", "--processor", help="Number of processors allowed (Optional)", default=os.cpu_count())
    parser.add_argument("-w", "--working_dir", help="Ou(Optional)", default=None)
    parser.add_argument("-o", "--output_dir", help="Output directory (Optional)", default=None)
    parser.add_argument("-d", "--debug", help="Ou (Optional)", action="store_true")
    parser.add_argument("-v", "--verbose", help="Verbose mode, will show interactive progress bar (Optional)",
                        action="store_true")
    # parser.add_argument('-c', '--checkpoint', help="Ou (Optional)", action="store_true")
    parser.add_argument("-c", "--checkpoint", help="Ou (Optional)", nargs='?', const="auto", default=None)
    parser.add_argument("-fc", "--force", help="Ou (Optional)", action="store_true")
    # TODO write help information

    parsed_args = parser.parse_args()

    genome_file = parsed_args.genome_file
    genome_name = parsed_args.genome_name
    species = parsed_args.species
    GRF_path = parsed_args.grfp
    GRF_mode = parsed_args.mode
    checkpoint_input = parsed_args.checkpoint

    TIR_length = int(parsed_args.length)
    cpu_cores = int(parsed_args.processor)
    working_dir = parsed_args.working_dir
    output_dir = parsed_args.output_dir
    if output_dir is None:
        output_dir = os.path.dirname(genome_file)

    flag_debug = parsed_args.debug
    flag_verbose = parsed_args.verbose
    # flag_checkpoint = parsed_args.checkpoint
    flag_force = parsed_args.force

    # Transforming the possible relative path into absolute path
    genome_file = os.path.abspath(genome_file)
    GRF_path = os.path.abspath(GRF_path)
    output_dir = os.path.abspath(output_dir)

    if checkpoint_input is not None and checkpoint_input != "auto":
        checkpoint_input = os.path.abspath(checkpoint_input)
    # ==================================================================================================================

    TIRLearner_instance = TIRLearner(genome_file, genome_name, species, TIR_length,
                                     working_dir, output_dir, GRF_path, cpu_cores, GRF_mode,
                                     checkpoint_input, flag_verbose, flag_debug, flag_force)
