#!/usr/bin/env python3
# Tianyu Lu (tlu83@wisc.edu)
# 2023-08-09

import os
import tempfile
import shutil
import argparse

import M1_1_blast_reference
import M1_2_full_coverage
import M2_3A_blast_reference
import M2_3B_eighty_similarity

import run_GRF
import process_GRFmite
import prepare_data
import CNN_predict
import get_fasta_sequence
import check_TIR_TSD
import post_processing


def execute_M1(args):
    print("############################################################ Module 1 Begin "
          "###########################################################")

    module = "Module1"
    # # Create genome_name directory if it doesn't exist
    # os.makedirs(os.path.join(dir, module), exist_ok=True)
    # os.chdir(os.path.join(dir, module))

    # Module 1, Step 1: Blast Genome against Reference Library
    M1_1_blast_reference.execute(args)

    # Module 1, Step 2: Select 100% coverage entries from Blast results
    df_full_cov = M1_2_full_coverage.execute(args)

    # Module 1, Step 3: Making blastDB and get candidate FASTA sequences
    print("Module 1, Step 3: Making blastDB and get candidate FASTA sequences")
    df_fasta = get_fasta_sequence.execute(args, df_full_cov)
    del df_full_cov

    # Module 1, Step 4: Check TIR and TSD
    print("Module 1, Step 4: Check TIR and TSD")
    df_TIR_TSD = check_TIR_TSD.execute(df_fasta, module)
    del df_fasta

    print("############################################################ Module 1 Finished "
          "########################################################")
    return df_TIR_TSD


def execute_M2(args):
    print("############################################################ Module 2 Begin "
          "###########################################################")

    module = "Module2"
    # os.makedirs(os.path.join(dir, module), exist_ok=True)
    # os.chdir(os.path.join(dir, module))

    # Module 2, Step 1: Split Genome and Run GRF program to find Inverted Repeats
    print("Module 2, Step 1: Run GRF program to find Inverted Repeats")
    run_GRF.execute(args)

    # Module 2, Step 2: Process GRF results
    print("Module 2, Step 2: Process GRF results")
    process_GRFmite.execute(args)

    # Module 2, Step 3A: GRF result blast reference sequences
    M2_3A_blast_reference.execute(args)

    # Module 2, Step 3B: Select 80% similar entries from blast results
    df_80_sim = M2_3B_eighty_similarity.execute(args)

    # Module 2, Step 4: Get FASTA sequences from 80% similarity
    print("Module 2, Step 4: Get FASTA sequences from 80% similarity")
    df_fasta = get_fasta_sequence.execute(args, df_80_sim)
    del df_80_sim

    # Module 2, Step 5: Check TIR and TSD
    print("Module 2, Step 5: Check TIR and TSD")
    df_TIR_TSD = check_TIR_TSD.execute(df_fasta, module)

    print("############################################################ Module 2 Finished "
          "########################################################")
    return df_TIR_TSD, df_fasta


def execute_M3(args, df_homo):
    print("############################################################ Module 3 Begin "
          "###########################################################")

    module = "Module3"
    # os.makedirs(os.path.join(dir, module), exist_ok=True)
    # os.chdir(os.path.join(dir, module))

    # Module 3, Step 1: Prepare Data
    print("Module 3, Step 1: Prepare Data")
    df_non_homo = prepare_data.execute(args, df_homo)

    # Module 3, Step 2: ML (CNN) prediction
    print("Module 3, Step 2: ML (CNN) prediction")
    df_pred = CNN_predict.execute(args, df_non_homo)
    del df_non_homo

    # Module 3, Step 3: Get FASTA sequences from ML prediction
    print("Module 3, Step 3: Get FASTA sequences from ML prediction")
    df_fasta = get_fasta_sequence.execute(args, df_pred)

    # Module 3, Step 4: Check TIR and TSD
    print("Module 3, Step 4: Check TIR and TSD")
    df_TIR_TSD = check_TIR_TSD.execute(df_fasta, module)
    del df_fasta

    print("############################################################ Module 3 Finished "
          "########################################################")
    return df_TIR_TSD


def execute_M3N(args):
    print("########################################################## Module 3 New Begin "
          "#########################################################")

    module = "Module3_New"
    # os.makedirs(os.path.join(dir, module), exist_ok=True)
    # os.chdir(os.path.join(dir, module))

    # Module 3N, Step 1: Split Genome and Run GRF program to find Inverted Repeats
    print("Module 3N, Step 1: Run GRF program to find Inverted Repeats")
    run_GRF.execute(args)

    # Module 3N, Step 2: Process GRF results
    print("Module 3N, Step 2: Process GRF results")
    process_GRFmite.execute(args)

    # Module 3N, Step 3: Prepare Data
    print("Module 3N, Step 3: Prepare Data")
    df_all = prepare_data.execute(args)

    # Module 3N, Step 4: ML (CNN) prediction
    print("Module 3N, Step 4: ML (CNN) prediction")
    df_pred = CNN_predict.execute(args, df_all)
    del df_all

    # Module 3N, Step 5: Get FASTA sequences from ML prediction
    print("Module 3N, Step 5: Get FASTA sequences from ML prediction")
    df_fasta = get_fasta_sequence.execute(args, df_pred)
    del df_pred

    # Module 3N, Step 6: Check TIR and TSD
    print("Module 3N, Step 6: Check TIR and TSD")
    df_TIR_TSD = check_TIR_TSD.execute(df_fasta, module)
    del df_fasta

    print("########################################################## Module 3 New Finished "
          "######################################################")
    return df_TIR_TSD


if __name__ == "__main__":
    # ================================================ argument parsing ================================================
    parser = argparse.ArgumentParser(prog="TIR-Learner")
    parser.add_argument("-f", "--genome_file", help="Genome file in fasta format", required=True)
    parser.add_argument("-n", "--genome_name", help="Genome name", required=True)
    parser.add_argument("-s", "--species", help="One of the following: \"Maize\", \"Rice\" or \"Others\"",
                        required=True)
    parser.add_argument("-c", "--CNN_path", help="Path to Tensorflow CNN SavedModel "
                                                 "(Do not include the SavedModel folder name)", required=True)
    parser.add_argument("-g", "--GRF_path", help="Path to GRF program", required=True)
    parser.add_argument("-l", "--TIR_length", help="Max length of TIR (Optional)", default=5000)
    parser.add_argument("-t", "--processor", help="Number of processor (Optional)", default=os.cpu_count())
    parser.add_argument("-o", "--output_dir", help="Output directory (Optional)", default=None)
    parser.add_argument('-d', '--debug', help="Output each module's result in csv file (Optional)", action="store_true")
    parsed_args = parser.parse_args()

    genome_file = parsed_args.genome_file
    genome_name = parsed_args.genome_name
    species = parsed_args.species
    CNN_path = parsed_args.CNN_path
    grfp = parsed_args.GRF_path

    length = parsed_args.TIR_length
    t = parsed_args.processor
    output_dir = parsed_args.output_dir
    flag_debug = parsed_args.debug
    # ==================================================================================================================

    temp_dir = tempfile.mkdtemp()
    genome_file_hard_link = os.path.join(temp_dir, "genome_file_hard_link.fa.lnk")
    os.link(genome_file, genome_file_hard_link)
    os.chdir(temp_dir)

    args_list = [genome_file_hard_link, genome_name, CNN_path, t, species, grfp, length]
    # index:                0                 1         2      3     4      5      6

    df_list = None
    if species == "Rice" or species == "Maize":
        df_M1 = execute_M1(args_list)
        df_M2, df_homo = execute_M2(args_list)
        df_M3 = execute_M3(args_list, df_homo)
        df_list = [df_M1, df_M2, df_M3]
    else:
        df_M3N = execute_M3N(args_list)
        df_list = [df_M3N]

    if output_dir is None or output_dir == "None":
        output_dir = os.path.dirname(genome_file)

    if flag_debug:
        for i, df in enumerate(df_list, start=1):
            df.to_csv(os.path.join(output_dir, f"debug_{i}.csv"), sep='\t')
    post_processing.execute(args_list, df_list, output_dir)

    # subprocess.Popen(["rm", "-f", f"{genome_file}.fai"])
    # subprocess.Popen(["rm", "-rf", temp_dir])
    shutil.rmtree(temp_dir)
