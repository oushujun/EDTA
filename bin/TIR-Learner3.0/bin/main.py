#!/usr/bin/env python3
# Tianyu Lu (tlu83@wisc.edu)
# 2023-09-22
import datetime
import json
import os
import re
import shutil
import tempfile

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

import M1_1_blast_reference
import M1_2_full_coverage
import M2_3_blast_reference
import M2_4_eighty_similarity

import prog_const
import run_GRF
import process_GRFmite
import prepare_data
import CNN_predict
import get_fasta_sequence
import check_TIR_TSD
import post_processing


class TIRLearner:
    def __init__(self, genome_file: str, genome_name: str, species: str, TIR_length: int,
                 working_dir: str, output_dir: str, GRF_path: str, cpu_cores: int, GRF_mode: str,
                 flag_verbose: bool, flag_debug: bool, flag_checkpoint: bool, flag_force: bool):

        self.genome_file = genome_file
        self.genome_name = genome_name
        self.working_dir = working_dir
        self.output_dir = output_dir
        self.species = species
        self.TIR_length = TIR_length
        self.GRF_path = GRF_path
        self.cpu_cores = cpu_cores
        self.GRF_mode = GRF_mode
        self.flag_verbose = flag_verbose
        self.flag_debug = flag_debug
        self.flag_checkpoint = flag_checkpoint
        self.flag_force = flag_force

        self.processedGRFmite_file = f"{self.genome_name}{prog_const.spliter}processedGRFmite.fa"
        if self.flag_debug or self.flag_checkpoint:
            timestamp_now_iso8601 = datetime.datetime.now().strftime("%Y-%m-%dT%H-%M-%S")
            self.checkpoint_folder = os.path.join(self.output_dir, f"TIR-Learner_v3_checkpoint_{timestamp_now_iso8601}")
            os.makedirs(self.checkpoint_folder)

        self.df_list = None
        self.genome_file_stat = {"file_size_gib": -0.1, "num": -1,
                                 "short_seq_num": 0, "short_seq_perc": -0.1,
                                 "total_len": 0, "avg_len": -1}
        self.current_step = [0, 0]
        self.working_df_dict = {}

        self.execute()

    def execute(self):
        self.mount_working_dir()
        self.load_checkpoint_file()
        self.check_fasta_file()
        # print(os.getcwd())  # TODO ONLY FOR DEBUG REMOVE AFTER FINISHED

        if self.species == "rice" or self.species == "maize":
            self.execute_M1()
            self.execute_M2()
            self.execute_M3()
            self.df_list = [self.working_df_dict["m1"], self.working_df_dict["m2"], self.working_df_dict["m3"]]
        else:
            self.execute_M4()
            self.df_list = [self.working_df_dict["m4"]]

        post_processing.execute(self)
        shutil.rmtree(self.working_dir)

    def check_fasta_file(self):
        # names = [record.id for record in SeqIO.parse(self.genome_file, "fasta")]
        print("Doing pre-scan for genome file...")
        self.genome_file_stat["file_size_gib"] = os.path.getsize(self.genome_file) / (2 ** 10) ** 3
        # if file_size_gib > 1.0:
        #     print(f"File size {file_size_gib} GiB > 1 GiB, pre-scan might be slow.")
        try:
            self.genome_file_stat["num"] = len(list(SeqIO.index(self.genome_file, "fasta")))
        except Exception:
            raise SystemExit(f"Duplicate sequence name occurs in the fasta file.")

        records = []
        for record in SeqIO.parse(self.genome_file, "fasta"):
            record.seq = record.seq.upper()
            sequence_str = str(record.seq)

            if set(sequence_str) > {'A', 'C', 'G', 'T', 'N'}:
                print(f"Unknown character exist in sequence {record.id} will be replaced by \'N\'")
                record.seq = Seq(re.sub("[^ACGTN]", "N", sequence_str))

            if prog_const.spliter in record.id:
                print((f"Sequence name \"{record.id}\" has reserved string \"{prog_const.spliter}\", "
                       "which makes it incompatible with TIR-Learner and will be replaced with \'_\'."))
                record.id = record.id.replace(prog_const.spliter, "_")

            if len(sequence_str) < prog_const.short_seq_len:
                self.genome_file_stat["short_seq_num"] += 1

            self.genome_file_stat["total_len"] += len(sequence_str)
            records.append(record)
        checked_genome_file = f"{self.genome_name}_checked.fa"
        SeqIO.write(records, checked_genome_file, "fasta")
        self.genome_file_stat["short_seq_perc"] = (self.genome_file_stat["short_seq_num"] /
                                                   self.genome_file_stat["num"])
        self.genome_file_stat["avg_len"] = self.genome_file_stat["total_len"] // self.genome_file_stat["num"]

        self.GRF_execution_mode_check()
        print("Genome file scan finished!")
        print(f"  File name: {os.path.basename(self.genome_file)}")
        print(f"  File size: {round(self.genome_file_stat['file_size_gib'], 4)} GiB")
        print(f"  Number of sequences: {self.genome_file_stat['num']}")
        print(f"  Number of short sequences: {self.genome_file_stat['short_seq_num']}")
        print(f"  Percentage of short sequences: {self.genome_file_stat['short_seq_perc']}")
        print(f"  Average sequence length: {self.genome_file_stat['avg_len'] // 1000} k")
        self.genome_file = os.path.abspath(checked_genome_file)

    def mount_working_dir(self):
        if self.working_dir is None:
            self.working_dir = tempfile.mkdtemp()
        # self.load_genome_file()
        self.working_dir = os.path.abspath(self.working_dir)
        os.chdir(self.working_dir)

    # def load_genome_file(self):
    #     genome_file_soft_link = os.path.join(self.execution_dir, "genome_file_soft_link.fa.lnk")
    #     os.symlink(self.genome_file, genome_file_soft_link)
    #     self.genome_file = genome_file_soft_link

    def get_newest_checkpoint_folder(self, search_dir: str):
        checkpoint_folders = [f for f in os.listdir(search_dir) if f.startswith("TIR-Learner_v3_checkpoint_")]
        try:
            checkpoint_folders.remove(self.checkpoint_folder)
        except ValueError:
            pass
        checkpoint_folders = sorted(checkpoint_folders)

        if len(checkpoint_folders) == 0:
            return None
        return os.path.join(self.output_dir, checkpoint_folders[-1])

    def load_checkpoint_file(self):
        if not self.flag_checkpoint:
            return

        # Search both the output directory and the genome file directory
        checkpoint_folder = self.get_newest_checkpoint_folder(self.output_dir)
        if checkpoint_folder is None:
            genome_file_directory = os.path.dirname(self.genome_file)
            checkpoint_folder = self.get_newest_checkpoint_folder(genome_file_directory)
        if checkpoint_folder is None:
            print("Unable to find checkpoint file. Will skip loading checkpoint and start from the very beginning.")
            # self.flag_checkpoint = False
            return

        checkpoint_info_file = os.path.join(checkpoint_folder, "info.txt")
        # print(checkpoint_info_file) # TODO only for debug
        if not os.path.exists(checkpoint_info_file) or os.path.getsize(checkpoint_info_file) == 0:
            print("Checkpoint file invalid. Will skip loading checkpoint and start from the very beginning.")
            # self.flag_checkpoint = False
            return

        with open(checkpoint_info_file) as f:
            timestamp_iso8601 = f.readline().rstrip()

            self.current_step = list(map(int, f.readline().rstrip().split(",")))
            module = self.current_step[0]
            step = self.current_step[1]

            # df_current_file = f.readline().rstrip()
            working_df_filename_dict = json.loads(f.readline().rstrip())
            for k, v in working_df_filename_dict.items():
                self.working_df_dict[k] = pd.read_csv(os.path.join(checkpoint_folder, v),
                                                      sep='\t', header=0, engine='c', memory_map=True)
            # self.df_current = pd.read_csv(df_current_file, sep='\t', header=0, engine='c', memory_map=True)

            if step == 2 and module in (2, 4):
                # Load processedGRFmite checkpoint file for module 2 step 2 or module 4 step 2
                processedGRFmite_checkpoint_file = os.path.join(checkpoint_folder, self.processedGRFmite_file)
                if os.path.exists(processedGRFmite_checkpoint_file):
                    shutil.copy(processedGRFmite_checkpoint_file, self.processedGRFmite_file)

            print(("Successfully loaded checkpoint:\n"
                   f"  Time: {timestamp_iso8601}\n"
                   f"  Module: {module}\n"
                   f"  Step: {step}"))

    def save_checkpoint_file(self):
        if not (self.flag_debug or self.flag_checkpoint):
            return

        # print(self.current_step) # TODO debug only
        module = self.current_step[0]
        step = self.current_step[1]
        timestamp_now_iso8601 = datetime.datetime.now().strftime("%Y-%m-%dT%H-%M-%S")
        # checkpoint_file_name = f"module_{module}_step_{step}_{timestamp_now_iso8601}.csv"
        working_df_filename_dict = {k: f"{k}_module_{module}_step_{step}_{timestamp_now_iso8601}.csv"
                                    for k in self.working_df_dict.keys()}
        # print(working_df_filename_dict) # TODO debug only
        for k, v in self.working_df_dict.items():
            v.to_csv(os.path.join(self.checkpoint_folder, working_df_filename_dict[k]),
                     index=False, header=True, sep="\t")
        # self.df_current.to_csv(os.path.join(self.checkpoint_folder, checkpoint_file_name),
        #                        index=False, header=False, sep="\t")
        with open(os.path.join(self.checkpoint_folder, "info.txt"), "w") as f:
            f.write(timestamp_now_iso8601 + '\n')
            f.write(f"{module},{step}\n")
            # f.write(checkpoint_file_name)
            f.write(json.dumps(working_df_filename_dict))
            f.write('\n')

        # print(os.listdir(self.checkpoint_folder)) # TODO debug only
        if not self.flag_debug:
            # shutil.rmtree(self.checkpoint_folder)
            remove_file_set = (set(os.listdir(self.checkpoint_folder)) -
                               set(working_df_filename_dict.values()) -
                               {"info.txt"})
            # print(remove_file_set) # TODO debug only
            for f in remove_file_set:
                os.remove(os.path.join(self.checkpoint_folder, f))

    def module_step_execution_check(self, executing_module: int, executing_step: int) -> bool:
        return (not self.flag_checkpoint or self.current_step[0] < executing_module or
                (self.current_step[0] == executing_module and self.current_step[1] < executing_step))

    def GRF_execution_mode_check(self):
        if self.flag_force or self.GRF_mode in ("native", "boost"):
            return

        if self.genome_file_stat["num"] <= prog_const.general_split_num_threshold:
            if self.GRF_mode == "smart":
                print("  \"native\" mode is selected due to insufficient number of sequences.")
            else:
                print(f"   Number of sequences insufficient "
                      f"(expect >= {prog_const.general_split_num_threshold}"
                      f" but actually got {self.genome_file_stat['num']}), "
                      f"{self.GRF_mode} mode unneeded, redirect to \"native\" mode.")
            self.GRF_mode = "native"
            return

        # "mix" mode or "smart" mode
        drop_seq_len = int(self.TIR_length) + 500
        if drop_seq_len >= prog_const.short_seq_len:
            if self.GRF_mode == "mix":
                print("  Short sequence does not exist after dropping, "
                      "\"mix\" mode unneeded, redirect to \"native\" mode.")
            else:
                print("  \"native\" mode is selected due to short sequence does not exist after dropping.")
            self.GRF_mode = "native"
            return

        if self.genome_file_stat["short_seq_perc"] < prog_const.mix_split_percent_threshold:
            if self.GRF_mode == "mix":
                print(f"  Percentage of short sequences insufficient "
                      f"(expect >= {prog_const.mix_split_percent_threshold * 100}%"
                      f" but actually got {self.genome_file_stat['short_seq_perc'] * 100}%), "
                      f"\"mix\" mode unneeded, redirect to \"native\" mode")
            else:
                print("  \"native\" mode is selected due to low percentage of short sequences.")
            self.GRF_mode = "native"
            return

        if self.genome_file_stat["short_seq_perc"] > 1 - prog_const.mix_split_percent_threshold:
            if self.GRF_mode == "mix":
                print(f"  Percentage of short sequences too high "
                      f"(expect < {(1 - prog_const.mix_split_percent_threshold) * 100}%"
                      f" but actually got {self.genome_file_stat['short_seq_perc'] * 100}%), "
                      f"\"mix\" mode inappropriate, redirect to \"boost\" mode")
            else:
                print("  \"boost\" mode is selected due to high percentage of short sequences.")
            self.GRF_mode = "boost"
            return

        if self.GRF_mode == "smart":
            print("  \"mix\" mode is selected.")
            self.GRF_mode = "mix"

    def execute_M1(self):
        print("############################################################ Module 1 Begin "
              "###########################################################")

        module = "Module1"
        # # Create genome_name directory if it doesn't exist
        # os.makedirs(os.path.join(dir, module), exist_ok=True)
        # os.chdir(os.path.join(dir, module))

        # Module 1, Step 1: Blast Genome against Reference Library
        if self.module_step_execution_check(1, 1):
            M1_1_blast_reference.execute(self)
            self.current_step = [1, 1]

        # Module 1, Step 2: Select 100% coverage entries from Blast results
        if self.module_step_execution_check(1, 2):
            self.working_df_dict["base"] = M1_2_full_coverage.execute(self)
            self.current_step = [1, 2]
            self.save_checkpoint_file()

        # Module 1, Step 3: Making blastDB and get candidate FASTA sequences
        if self.module_step_execution_check(1, 3):
            print("Module 1, Step 3: Making blastDB and get candidate FASTA sequences")
            self.working_df_dict["base"] = get_fasta_sequence.execute(self)
            self.current_step = [1, 3]
            self.save_checkpoint_file()

        # Module 1, Step 4: Check TIR and TSD
        if self.module_step_execution_check(1, 4):
            print("Module 1, Step 4: Check TIR and TSD")
            self.working_df_dict["base"] = check_TIR_TSD.execute(self, module)
            self.current_step = [1, 4]
            self.save_checkpoint_file()

        # Module 1, Step 5: Save module result
        if self.module_step_execution_check(1, 5):
            self.working_df_dict["m1"] = self.working_df_dict["base"]
            del self.working_df_dict["base"]
            self.current_step = [1, 5]
            self.save_checkpoint_file()

        print("############################################################ Module 1 Finished "
              "########################################################")
        # return self.df_current.copy()

    def execute_M2(self):
        print("############################################################ Module 2 Begin "
              "###########################################################")

        module = "Module2"
        # os.makedirs(os.path.join(dir, module), exist_ok=True)
        # os.chdir(os.path.join(dir, module))

        # Module 2, Step 1: Split Genome and Run GRF program to find Inverted Repeats
        if self.module_step_execution_check(2, 1):
            # Checkpoint saving for this step is currently not available
            print("Module 2, Step 1: Run GRF program to find Inverted Repeats")
            run_GRF.execute(self)
            # self.current_step = [2, 1]

        # Module 2, Step 2: Process GRF results
        print("Module 2, Step 2: Process GRF results")
        if self.module_step_execution_check(2, 2):
            process_GRFmite.execute(self)
            self.current_step = [2, 2]
            if self.flag_debug or self.flag_checkpoint:
                shutil.copy(self.processedGRFmite_file,
                            os.path.join(self.checkpoint_folder, self.processedGRFmite_file))

        # Module 2, Step 3: GRF result blast reference sequences
        if self.module_step_execution_check(2, 3):
            # Checkpoint saving for this step is currently not available
            M2_3_blast_reference.execute(self)
            self.current_step = [2, 3]

        # Module 2, Step 4: Select 80% similar entries from blast results
        if self.module_step_execution_check(2, 4):
            self.working_df_dict["base"] = M2_4_eighty_similarity.execute(self)
            self.current_step = [2, 4]
            self.save_checkpoint_file()

        # Module 2, Step 5: Get FASTA sequences from 80% similarity
        if self.module_step_execution_check(2, 5):
            print("Module 2, Step 5: Get FASTA sequences from 80% similarity")
            self.working_df_dict["base"] = get_fasta_sequence.execute(self)
            self.working_df_dict["m2_homo"] = self.working_df_dict["base"].copy()
            self.current_step = [2, 5]
            self.save_checkpoint_file()

        # Module 2, Step 6: Check TIR and TSD
        if self.module_step_execution_check(2, 6):
            print("Module 2, Step 6: Check TIR and TSD")
            self.working_df_dict["base"] = check_TIR_TSD.execute(self, module)
            self.current_step = [2, 6]
            self.save_checkpoint_file()

        # Module 2, Step 7: Save module result
        if self.module_step_execution_check(2, 7):
            self.working_df_dict["m2"] = self.working_df_dict["base"]
            del self.working_df_dict["base"]
            self.current_step = [2, 7]
            self.save_checkpoint_file()

        print("############################################################ Module 2 Finished "
              "########################################################")
        # return self.df_current.copy(), df_fasta

    def execute_M3(self):
        print("############################################################ Module 3 Begin "
              "###########################################################")

        module = "Module3"
        # os.makedirs(os.path.join(dir, module), exist_ok=True)
        # os.chdir(os.path.join(dir, module))

        # Module 3, Step 1: Prepare Data
        if self.module_step_execution_check(3, 1):
            print("Module 3, Step 1: Prepare Data")
            self.working_df_dict["base"] = prepare_data.execute(self, self.working_df_dict["m2_homo"])
            del self.working_df_dict["m2_homo"]
            self.current_step = [3, 1]
            self.save_checkpoint_file()

        # Module 3, Step 2: ML (CNN) prediction
        if self.module_step_execution_check(3, 2):
            print("Module 3, Step 2: ML (CNN) prediction")
            self.working_df_dict["base"] = CNN_predict.execute(self)
            self.current_step = [3, 2]
            self.save_checkpoint_file()

        # Module 3, Step 3: Get FASTA sequences from ML prediction
        if self.module_step_execution_check(3, 3):
            print("Module 3, Step 3: Get FASTA sequences from ML prediction")
            self.working_df_dict["base"] = get_fasta_sequence.execute(self)
            self.current_step = [3, 3]
            self.save_checkpoint_file()

        # Module 3, Step 4: Check TIR and TSD
        if self.module_step_execution_check(3, 4):
            print("Module 3, Step 4: Check TIR and TSD")
            self.working_df_dict["base"] = check_TIR_TSD.execute(self, module)
            self.current_step = [3, 4]
            self.save_checkpoint_file()

        # Module 3, Step 5: Save module result
        if self.module_step_execution_check(3, 5):
            self.working_df_dict["m3"] = self.working_df_dict["base"]
            del self.working_df_dict["base"]
            self.current_step = [3, 5]
            self.save_checkpoint_file()

        print("############################################################ Module 3 Finished "
              "########################################################")
        # return df_TIR_TSD

    def execute_M4(self):
        print("########################################################## Module 4 Begin "
              "#########################################################")

        module = "Module4"
        # os.makedirs(os.path.join(dir, module), exist_ok=True)
        # os.chdir(os.path.join(dir, module))

        # Module 4, Step 1: Split Genome and Run GRF program to find Inverted Repeats
        if self.module_step_execution_check(4, 1):
            print("Module 4, Step 1: Run GRF program to find Inverted Repeats")
            # Checkpoint saving for this step is currently not available
            run_GRF.execute(self)
            # self.current_step = [4, 1]

        # Module 4, Step 2: Process GRF results
        if self.module_step_execution_check(4, 2):
            print("Module 4, Step 2: Process GRF results")
            process_GRFmite.execute(self)
            self.current_step = [4, 2]
            if self.flag_debug or self.flag_checkpoint:
                shutil.copy(self.processedGRFmite_file,
                            os.path.join(self.checkpoint_folder, self.processedGRFmite_file))

        # Module 4, Step 3: Prepare Data
        if self.module_step_execution_check(4, 3):
            print("Module 4, Step 3: Prepare Data")
            self.working_df_dict["base"] = prepare_data.execute(self)
            self.current_step = [4, 3]
            self.save_checkpoint_file()

        # Module 4, Step 4: ML (CNN) prediction
        print("Module 4, Step 4: ML (CNN) prediction")
        if self.module_step_execution_check(4, 4):
            self.working_df_dict["base"] = CNN_predict.execute(self)
            self.current_step = [4, 4]
            self.save_checkpoint_file()

        # Module 4, Step 5: Get FASTA sequences from ML prediction
        if self.module_step_execution_check(4, 5):
            print("Module 4, Step 5: Get FASTA sequences from ML prediction")
            self.working_df_dict["base"] = get_fasta_sequence.execute(self)
            self.current_step = [4, 5]
            self.save_checkpoint_file()

        # Module 4, Step 6: Check TIR and TSD
        if self.module_step_execution_check(4, 6):
            print("Module 4, Step 6: Check TIR and TSD")
            self.working_df_dict["base"] = check_TIR_TSD.execute(self, module)
            self.current_step = [4, 6]
            self.save_checkpoint_file()

        # Module 4, Step 7: Save module result
        if self.module_step_execution_check(4, 7):
            self.working_df_dict["m4"] = self.working_df_dict["base"]
            del self.working_df_dict["base"]
            self.current_step = [4, 7]
            self.save_checkpoint_file()

        print("########################################################## Module 4 Finished "
              "######################################################")
        # return self.df_current.copy()

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
