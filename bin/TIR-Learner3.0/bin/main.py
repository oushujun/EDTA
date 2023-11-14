#!/usr/bin/env python3
# Tianyu Lu (tlu83@wisc.edu)
# 2023-09-22

import os

import pandas as pd
from Bio import SeqIO

import M1_1_blast_reference
import M1_2_full_coverage
import M2_3A_blast_reference
import M2_3B_eighty_similarity

import prog_const
import run_GRF
import process_GRFmite
import prepare_data
import CNN_predict
import get_fasta_sequence
import check_TIR_TSD
import post_processing


class TIRLearner:
    def __init__(self, genome_file: str, genome_name: str, output_dir: str, species: str, TIR_length: int,
                 GRF_path: str, cpu_cores: int, GRF_mode: str, flag_verbose: bool, flag_debug: bool):

        self.genome_file = genome_file
        self.genome_name = genome_name
        self.output_dir = output_dir
        self.species = species
        self.TIR_length = TIR_length
        self.GRF_path = GRF_path
        self.cpu_cores = cpu_cores
        self.GRF_mode = GRF_mode
        self.flag_verbose = flag_verbose
        self.flag_debug = flag_debug

        self.processedGRFmite_file = f"{self.genome_name}{prog_const.spliter}processedGRFmite.fa"

        self.df_list = None
        self.genome_num = None

        self.execute()

    def execute(self):
        self.check_fasta_file()
        print(os.getcwd())

        if self.species == "rice" or self.species == "maize":
            df_M1 = self.execute_M1()
            df_M2, df_homo = self.execute_M2()
            df_M3 = self.execute_M3(df_homo)
            self.df_list = [df_M1, df_M2, df_M3]
        else:
            df_M3N = self.execute_M3N()
            self.df_list = [df_M3N]

        self.debug()
        post_processing.execute(self)

    def check_fasta_file(self):
        # names = [record.id for record in SeqIO.parse(self.genome_file, "fasta")]
        names = list(SeqIO.index(self.genome_file, "fasta"))

        error_name = next((name for name in names if prog_const.spliter in name), None)
        # if any(self.spliter in name for name in names):
        if error_name is not None:
            raise SystemExit((f"Sequence name \"{error_name}\" has reserved string \"{prog_const.spliter}\", "
                              "which makes it incompatible with TIR-Learner, "
                              "please remove this reserved string from the name."))

        error_name = next((name for name in names if names.count(name) > 1), None)
        if error_name is not None:
            raise SystemExit(f"Duplicate sequence name \"{error_name}\" occurs in the fasta file.")

        self.genome_num = len(names)
        print(f"{self.genome_num} sequences detected in the input genome file.")

    def debug(self):
        if self.flag_debug:
            for i, df in enumerate(self.df_list, start=1):
                df.to_csv(os.path.join(self.output_dir, f"debug_{i}.csv"), sep='\t')

    def execute_M1(self):
        print("############################################################ Module 1 Begin "
              "###########################################################")

        module = "Module1"
        # # Create genome_name directory if it doesn't exist
        # os.makedirs(os.path.join(dir, module), exist_ok=True)
        # os.chdir(os.path.join(dir, module))

        # Module 1, Step 1: Blast Genome against Reference Library
        M1_1_blast_reference.execute(self)

        # Module 1, Step 2: Select 100% coverage entries from Blast results
        df_full_cov = M1_2_full_coverage.execute(self)

        # Module 1, Step 3: Making blastDB and get candidate FASTA sequences
        print("Module 1, Step 3: Making blastDB and get candidate FASTA sequences")
        df_fasta = get_fasta_sequence.execute(self, df_full_cov)
        del df_full_cov

        # Module 1, Step 4: Check TIR and TSD
        print("Module 1, Step 4: Check TIR and TSD")
        df_TIR_TSD = check_TIR_TSD.execute(self, df_fasta, module)
        del df_fasta

        print("############################################################ Module 1 Finished "
              "########################################################")
        return df_TIR_TSD

    def execute_M2(self):
        print("############################################################ Module 2 Begin "
              "###########################################################")

        module = "Module2"
        # os.makedirs(os.path.join(dir, module), exist_ok=True)
        # os.chdir(os.path.join(dir, module))

        # Module 2, Step 1: Split Genome and Run GRF program to find Inverted Repeats
        print("Module 2, Step 1: Run GRF program to find Inverted Repeats")
        run_GRF.execute(self)

        # Module 2, Step 2: Process GRF results
        print("Module 2, Step 2: Process GRF results")
        process_GRFmite.execute(self)

        # Module 2, Step 3A: GRF result blast reference sequences
        M2_3A_blast_reference.execute(self)

        # Module 2, Step 3B: Select 80% similar entries from blast results
        df_80_sim = M2_3B_eighty_similarity.execute(self)

        # Module 2, Step 4: Get FASTA sequences from 80% similarity
        print("Module 2, Step 4: Get FASTA sequences from 80% similarity")
        df_fasta = get_fasta_sequence.execute(self, df_80_sim)
        del df_80_sim

        # Module 2, Step 5: Check TIR and TSD
        print("Module 2, Step 5: Check TIR and TSD")
        df_TIR_TSD = check_TIR_TSD.execute(self, df_fasta, module)

        print("############################################################ Module 2 Finished "
              "########################################################")
        return df_TIR_TSD, df_fasta

    def execute_M3(self, df_homo: pd.DataFrame):
        print("############################################################ Module 3 Begin "
              "###########################################################")

        module = "Module3"
        # os.makedirs(os.path.join(dir, module), exist_ok=True)
        # os.chdir(os.path.join(dir, module))

        # Module 3, Step 1: Prepare Data
        print("Module 3, Step 1: Prepare Data")
        df_non_homo = prepare_data.execute(self, df_homo)

        # Module 3, Step 2: ML (CNN) prediction
        print("Module 3, Step 2: ML (CNN) prediction")
        df_pred = CNN_predict.execute(self, df_non_homo)
        del df_non_homo

        # Module 3, Step 3: Get FASTA sequences from ML prediction
        print("Module 3, Step 3: Get FASTA sequences from ML prediction")
        df_fasta = get_fasta_sequence.execute(self, df_pred)

        # Module 3, Step 4: Check TIR and TSD
        print("Module 3, Step 4: Check TIR and TSD")
        df_TIR_TSD = check_TIR_TSD.execute(self, df_fasta, module)
        del df_fasta

        print("############################################################ Module 3 Finished "
              "########################################################")
        return df_TIR_TSD

    def execute_M3N(self):
        print("########################################################## Module 3 Begin "
              "#########################################################")

        module = "Module3"
        # os.makedirs(os.path.join(dir, module), exist_ok=True)
        # os.chdir(os.path.join(dir, module))

        # Module 3, Step 1: Split Genome and Run GRF program to find Inverted Repeats
        print("Module 3, Step 1: Run GRF program to find Inverted Repeats")
        run_GRF.execute(self)

        # Module 3, Step 2: Process GRF results
        print("Module 3, Step 2: Process GRF results")
        process_GRFmite.execute(self)

        # Module 3, Step 3: Prepare Data
        print("Module 3, Step 3: Prepare Data")
        df_all = prepare_data.execute(self)

        # Module 3, Step 4: ML (CNN) prediction
        print("Module 3, Step 4: ML (CNN) prediction")
        df_pred = CNN_predict.execute(self, df_all)
        del df_all

        # Module 3, Step 5: Get FASTA sequences from ML prediction
        print("Module 3, Step 5: Get FASTA sequences from ML prediction")
        df_fasta = get_fasta_sequence.execute(self, df_pred)
        del df_pred

        # Module 3, Step 6: Check TIR and TSD
        print("Module 3, Step 6: Check TIR and TSD")
        df_TIR_TSD = check_TIR_TSD.execute(self, df_fasta, module)
        del df_fasta

        print("########################################################## Module 3 Finished "
              "######################################################")
        return df_TIR_TSD


# def execute(args_list, output_dir, flag_debug):
#     genome_file = args_list[0]
#     species = args_list[4]
#
#     pre.check_sequence_names(genome_file)
#
#     if species == "rice" or species == "maize":
#         df_M1 = execute_M1(args_list)
#         df_M2, df_homo = execute_M2(args_list)
#         df_M3 = execute_M3(args_list, df_homo)
#         df_list = [df_M1, df_M2, df_M3]
#     else:
#         df_M3N = execute_M3N(args_list)
#         df_list = [df_M3N]
#
#     if flag_debug:
#         for i, df in enumerate(df_list, start=1):
#             df.to_csv(os.path.join(output_dir, f"debug_{i}.csv"), sep='\t')
#     post_processing.execute(args_list, df_list, output_dir)


# if __name__ == "__main__":
#     # ================================================ argument parsing ================================================
#     parser = argparse.ArgumentParser(prog="TIR-Learner")
#     parser.add_argument("-f", "--genome_file", help="Genome file in fasta format", required=True)
#     parser.add_argument("-n", "--genome_name", help="Genome name", required=True)
#     parser.add_argument("-s", "--species", help="One of the following: \"Maize\", \"Rice\" or \"Others\"",
#                         required=True)
#     parser.add_argument("-c", "--CNN_path", help="Path to Tensorflow CNN SavedModel "
#                                                  "(Do not include the SavedModel folder name)", required=True)
#     parser.add_argument("-g", "--GRF_path", help="Path to GRF program", required=True)
#     parser.add_argument("-l", "--TIR_length", help="Max length of TIR (Optional)", default=5000)
#     parser.add_argument("-t", "--processor", help="Number of processor (Optional)", default=os.cpu_count())
#     parser.add_argument("-o", "--output_dir", help="Output directory (Optional)", default=None)
#     parser.add_argument('-d', '--debug', help="Output each module's result in csv file (Optional)", action="store_true")
#     parsed_args = parser.parse_args()
#
#     genome_file = parsed_args.genome_file
#     genome_name = parsed_args.genome_name
#     species = parsed_args.species
#     CNN_path = parsed_args.CNN_path
#     grfp = parsed_args.GRF_path
#
#     length = parsed_args.TIR_length
#     t = parsed_args.processor
#     output_dir = parsed_args.output_dir
#     flag_debug = parsed_args.debug
#     # ==================================================================================================================
#
#     temp_dir = tempfile.mkdtemp()
#     genome_file_hard_link = os.path.join(temp_dir, "genome_file_hard_link.fa.lnk")
#     os.symlink(genome_file, genome_file_hard_link)
#     os.chdir(temp_dir)
#
#     args_list = [genome_file_hard_link, genome_name, CNN_path, t, species, grfp, length]
#     # index:                0                 1         2      3     4      5      6
#
#     df_list = None
#     if species == "Rice" or species == "Maize":
#         df_M1 = execute_M1(args_list)
#         df_M2, df_homo = execute_M2(args_list)
#         df_M3 = execute_M3(args_list, df_homo)
#         df_list = [df_M1, df_M2, df_M3]
#     else:
#         df_M3N = execute_M3N(args_list)
#         df_list = [df_M3N]
#
#     if output_dir is None or output_dir == "None":
#         output_dir = os.path.dirname(genome_file)
#
#     if flag_debug:
#         for i, df in enumerate(df_list, start=1):
#             df.to_csv(os.path.join(output_dir, f"debug_{i}.csv"), sep='\t')
#     post_processing.execute(args_list, df_list, output_dir)
#
#     # subprocess.Popen(["rm", "-f", f"{genome_file}.fai"])
#     # subprocess.Popen(["rm", "-rf", temp_dir])
#     shutil.rmtree(temp_dir)
