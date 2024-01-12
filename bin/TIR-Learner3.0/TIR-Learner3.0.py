#!/usr/bin/env python3
# Tianyu Lu (tlu83@wisc.edu)
# 2024-01-08

import argparse
import os
import shutil

import sys
sys.path.insert(0, f"{os.path.dirname(__file__)}/bin")

# Use if True to suppress the PEP8: E402 warning
if True:  # noqa: E402
    from bin.main import TIRLearner
    from bin import prog_const


if __name__ == "__main__":

    # ================================================ argument parsing ================================================
    parser = argparse.ArgumentParser(prog="TIR-Learner")
    parser.add_argument("-f", "--genome_file", help="Genome file in fasta format", required=True)
    parser.add_argument("-n", "--genome_name", help="Genome name (Optional)", default="TIR-Learner")
    parser.add_argument("-s", "--species", help="One of the following: \"maize\", \"rice\" or \"others\"",
                        required=True)
    parser.add_argument("-l", "--length", help="Max length of TIR (Optional)", default=5000)
    parser.add_argument("-t", "--processor", help="Number of processors allowed (Optional)", default=os.cpu_count())
    parser.add_argument('-m', '--mode', help=("Execution mode of GRF, one of the following: \"boost\", \"mix\""
                                              " or \"native\" (Optional)"), default="smart")
    parser.add_argument("-w", "--working_dir", help="Ou(Optional)", default=None)
    parser.add_argument("-o", "--output_dir", help="Output directory (Optional)", default=None)
    parser.add_argument("-c", "--checkpoint", help="Ou (Optional)", nargs='?', const="auto", default=None)
    parser.add_argument("-v", "--verbose", help="Verbose mode, will show interactive progress bar (Optional)",
                        action="store_true")
    parser.add_argument("-d", "--debug", help="Ou (Optional)", action="store_true")
    # parser.add_argument('-c', '--checkpoint', help="Ou (Optional)", action="store_true")
    parser.add_argument("--grf_path", help="Path to GRF program (Optional)",
                        default=os.path.dirname(shutil.which("grf-main")))
    parser.add_argument("--gt_path", help="Path to genometools program (Optional)",
                        default=os.path.dirname(shutil.which("gt")))
    parser.add_argument("-a", "--additional_args", help="Ou (Optional)", default="")
    # see prog_const for what additional args are acceptable
    # TODO write help information

    parsed_args = parser.parse_args()

    genome_file = parsed_args.genome_file
    genome_name = parsed_args.genome_name
    species = parsed_args.species

    TIR_length = int(parsed_args.length)
    cpu_cores = int(parsed_args.processor)
    GRF_mode = parsed_args.mode

    working_dir = parsed_args.working_dir
    output_dir = parsed_args.output_dir
    if output_dir is None:
        output_dir = os.path.dirname(genome_file)
    checkpoint_input = parsed_args.checkpoint

    flag_verbose = parsed_args.verbose
    flag_debug = parsed_args.debug

    GRF_path = parsed_args.grf_path.replace('"', "")
    gt_path = parsed_args.gt_path.replace('"', "")
    additional_args = prog_const.process_additional_args(parsed_args.additional_args.split(" "))
    if len(additional_args) != 0:
        print(f"INFO: Additional args: {additional_args} accepted.")

    # Transforming the possible relative path into absolute path
    genome_file = os.path.abspath(genome_file)
    output_dir = os.path.abspath(output_dir)
    GRF_path = os.path.abspath(GRF_path)
    gt_path = os.path.abspath(gt_path)

    if checkpoint_input is not None and checkpoint_input != "auto":
        checkpoint_input = os.path.abspath(checkpoint_input)
    # ==================================================================================================================

    TIRLearner_instance = TIRLearner(genome_file, genome_name, species, TIR_length,
                                     cpu_cores, GRF_mode, working_dir, output_dir, checkpoint_input,
                                     flag_verbose, flag_debug, GRF_path, gt_path, additional_args)
