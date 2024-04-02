#!/usr/bin/env python3
# Tianyu Lu (tlu83@wisc.edu)
# 2024-01-20

import os
import sys

sys.path.insert(0, f"{os.path.dirname(__file__)}/bin")

# Use if True to suppress the PEP8: E402 warning
if True:  # noqa: E402
    import argparse
    import shutil

    from bin.main import TIRLearner
    from bin.const import process_additional_args

if __name__ == "__main__":

    # ================================================ argument parsing ================================================
    parser = argparse.ArgumentParser(prog="TIR-Learner")
    parser.add_argument("-f", "--genome_file", help="Genome file in fasta format", required=True)
    parser.add_argument("-n", "--genome_name", help="Genome name (Optional)", default="TIR-Learner")
    parser.add_argument("-s", "--species", help="One of the following: \"maize\", \"rice\" or \"others\"",
                        required=True)
    parser.add_argument("-l", "--length", help="Max length of TIR (Optional)", default=5000)
    parser.add_argument("-t", "--processor", help="Number of processors allowed (Optional)", default=os.cpu_count())
    # TODO add pyboost, pystrict, gnu three parallel execution mode, also add more detailed help info
    parser.add_argument('-m', '--mode', help=("Parallel execution mode, one of the following: \"pyboost\", "
                                              "\"pystrict\" and \"gnup\" (Optional)"), default="pyboost")
    parser.add_argument("-w", "--working_dir", help="The path to the working directory (Optional). "
                                                    "An isolated sandbox directory for storing all the temporary files "
                                                    "will be created in the working directory. This sandbox directory "
                                                    "will only persist during the program execution. DO NOT TOUCH "
                                                    "THE SANDBOX DIRECTORY IF IT IS NOT FOR DEBUGGING!", default=None)
    parser.add_argument("-o", "--output_dir", help="Output directory (Optional)", default=None)
    parser.add_argument("-c", "--checkpoint_dir", help="The path to the checkpoint directory (Optional). "
                                                       "If not specified, the program will automatically search for it "
                                                       "in the genome file directory and the output directory.",
                        nargs='?', const="auto", default=None)
    parser.add_argument("-v", "--verbose", help="Verbose mode (Optional). "
                                                "Will show interactive progress bar and more execution details.",
                        action="store_true")
    parser.add_argument("-d", "--debug", help="Debug mode (Optional). If activated, data for all completed steps "
                                              "will be stored in the checkpoint file. Meanwhile, the temporary files "
                                              "in the working directory will also be kept.", action="store_true")
    parser.add_argument("--grf_path", help="Path to GRF program (Optional)",
                        default=os.path.dirname(shutil.which("grf-main")))
    parser.add_argument("--gt_path", help="Path to genometools program (Optional)",
                        default=os.path.dirname(shutil.which("gt")))
    parser.add_argument("-a", "--additional_args", help="Additional arguments (Optional). "
                                                        "See documentation for more details.", default="")
    # see prog_const for what additional args are acceptable

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
    checkpoint_input = parsed_args.checkpoint_dir

    flag_verbose = parsed_args.verbose
    flag_debug = parsed_args.debug

    GRF_path = parsed_args.grf_path.replace('"', "")
    gt_path = parsed_args.gt_path.replace('"', "")
    additional_args = process_additional_args(parsed_args.additional_args.split(" "))
    if len(additional_args) != 0:
        print(f"INFO: Additional args: {additional_args} captured.")

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
