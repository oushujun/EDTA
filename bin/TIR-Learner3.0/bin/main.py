#!/usr/bin/env python3
# Tianyu Lu (skyl@cs.wisc.edu)
# 2024-01-03
import datetime
import json
import os
import re
import shutil
import subprocess
import tempfile

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

import blast_reference
import process_homology
import run_TIRvish
import run_GRF
import process_de_novo_result
import prepare_data
import CNN_predict
import get_fasta_sequence
import check_TIR_TSD
import post_processing
import prog_const


def get_timestamp_now_utc_iso8601():
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H-%M-%SZ")


class TIRLearner:
    def __init__(self, genome_file: str, genome_name: str, species: str, TIR_length: int,
                 working_dir: str, output_dir: str, GRF_path: str, cpu_cores: int, GRF_mode: str,
                 checkpoint_input: str, flag_verbose: bool, flag_debug: bool, flag_force: bool):

        self.genome_file = genome_file
        self.genome_name = genome_name
        self.working_dir = working_dir
        self.output_dir = output_dir
        self.species = species
        self.TIR_length = TIR_length
        self.GRF_path = GRF_path
        self.cpu_cores = cpu_cores
        self.GRF_mode = GRF_mode
        self.checkpoint_input = checkpoint_input
        self.flag_verbose = flag_verbose
        self.flag_debug = flag_debug
        # self.flag_checkpoint = flag_checkpoint
        self.flag_force = flag_force

        self.processed_de_novo_result_file = f"{self.genome_name}{prog_const.spliter}processed_de_novo_result.fa"
        self.checkpoint_output = os.path.join(self.output_dir,
                                              f"TIR-Learner_v3_checkpoint_{get_timestamp_now_utc_iso8601()}")
        os.makedirs(self.checkpoint_output)

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
        self.pre_scan_fasta_file()
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
        if not self.flag_debug:
            shutil.rmtree(self.checkpoint_output)
        shutil.rmtree(self.working_dir)

    def pre_scan_fasta_file(self):
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
        print(f"  File size: {float('%.4g' % self.genome_file_stat['file_size_gib'])} GiB")
        print(f"  Number of sequences: {self.genome_file_stat['num']}")
        print(f"  Number of short sequences: {self.genome_file_stat['short_seq_num']}")
        print(f"  Percentage of short sequences: {self.genome_file_stat['short_seq_perc'] * 100} %")
        print(f"  Average sequence length: {self.genome_file_stat['avg_len'] // 1000} k")
        self.genome_file = os.path.abspath(checked_genome_file)

    def mount_working_dir(self):
        if self.working_dir is None:
            self.working_dir = tempfile.mkdtemp()
        # self.load_genome_file()
        else:
            os.makedirs(self.working_dir, exist_ok=True)
        self.working_dir = os.path.abspath(self.working_dir)
        os.chdir(self.working_dir)

    # def load_genome_file(self):
    #     genome_file_soft_link = os.path.join(self.execution_dir, "genome_file_soft_link.fa.lnk")
    #     os.symlink(self.genome_file, genome_file_soft_link)
    #     self.genome_file = genome_file_soft_link

    def get_newest_checkpoint_folder(self, search_dir: str) -> str:
        checkpoint_folders = [f for f in os.listdir(search_dir) if f.startswith("TIR-Learner_v3_checkpoint_")]
        # print(checkpoint_folders)  # TODO only for debug
        # print(self.checkpoint_folder)  # TODO only for debug
        try:
            checkpoint_folders.remove(os.path.basename(self.checkpoint_output))
        except ValueError:
            pass
        checkpoint_folders = sorted(checkpoint_folders)

        if len(checkpoint_folders) == 0:
            return "no_checkpoint_folder_found"
        return os.path.join(self.output_dir, checkpoint_folders[-1])

    def load_checkpoint_file(self):
        if self.checkpoint_input is None:
            return

        if self.checkpoint_input == "auto":
            # Search both the output directory and the genome file directory
            self.checkpoint_input = self.get_newest_checkpoint_folder(self.output_dir)
            if self.checkpoint_input == "no_checkpoint_folder_found":
                genome_file_directory = os.path.dirname(self.genome_file)
                self.checkpoint_input = self.get_newest_checkpoint_folder(genome_file_directory)
            if self.checkpoint_input == "no_checkpoint_folder_found":
                print("Unable to find checkpoint file. Will skip loading checkpoint and start from the very beginning.")
                # self.flag_checkpoint = False
                return

        # print(self.checkpoint_input)  # TODO only for debug
        checkpoint_info_file = os.path.join(self.checkpoint_input, "info.txt")
        # print(checkpoint_info_file)  # TODO only for debug
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
                self.working_df_dict[k] = pd.read_csv(os.path.join(self.checkpoint_input, v),
                                                      sep='\t', header=0, engine='c', memory_map=True)
            # self.df_current = pd.read_csv(df_current_file, sep='\t', header=0, engine='c', memory_map=True)

            # if step == 2 and module in (2, 4):
            #     # Load processedGRFmite checkpoint file for module 2 step 5 or module 4 step 2
            processed_de_novo_result_checkpoint_file = os.path.join(self.checkpoint_input,
                                                                    self.processed_de_novo_result_file)
            if os.path.exists(processed_de_novo_result_checkpoint_file):
                shutil.copy(processed_de_novo_result_checkpoint_file,
                            os.path.join(self.working_dir, self.processed_de_novo_result_file))
                shutil.copy(processed_de_novo_result_checkpoint_file,
                            os.path.join(self.checkpoint_output, self.processed_de_novo_result_file))

            print(("Successfully loaded checkpoint:\n"
                   f"  Time: {timestamp_iso8601}\n"
                   f"  Module: {module}\n"
                   f"  Step: {step}"))

    def save_checkpoint_file(self):
        # if not self.flag_debug and self.checkpoint_input is None:
        #     return

        # print(self.current_step) # TODO debug only
        module = self.current_step[0]
        step = self.current_step[1]
        # checkpoint_file_name = f"module_{module}_step_{step}_{timestamp_now_iso8601}.csv"
        working_df_filename_dict = {k: f"{k}_module_{module}_step_{step}_{get_timestamp_now_utc_iso8601()}.csv"
                                    for k in self.working_df_dict.keys()}
        # print(working_df_filename_dict) # TODO debug only
        for k, v in self.working_df_dict.items():
            v.to_csv(os.path.join(self.checkpoint_output, working_df_filename_dict[k]),
                     index=False, header=True, sep="\t")
        # self.df_current.to_csv(os.path.join(self.checkpoint_folder, checkpoint_file_name),
        #                        index=False, header=False, sep="\t")
        with open(os.path.join(self.checkpoint_output, "info.txt"), 'w') as f:
            f.write(get_timestamp_now_utc_iso8601() + '\n')
            f.write(f"{module},{step}\n")
            # f.write(checkpoint_file_name)
            f.write(json.dumps(working_df_filename_dict))
            f.write('\n')

        # print(os.listdir(self.checkpoint_folder)) # TODO debug only
        if not self.flag_debug:
            # shutil.rmtree(self.checkpoint_folder)
            remove_file_set = (set(os.listdir(self.checkpoint_output)) -
                               set(working_df_filename_dict.values()) -
                               {"info.txt", self.processed_de_novo_result_file})
            # print(remove_file_set) # TODO debug only
            for f in remove_file_set:
                # os.remove(os.path.join(self.checkpoint_output, f))
                subprocess.Popen(["unlink", os.path.join(self.checkpoint_output, f)])

    def save_processed_de_novo_result_checkpoint_file(self):
        shutil.copy(self.processed_de_novo_result_file,
                    os.path.join(self.checkpoint_output, self.processed_de_novo_result_file))
        with open(os.path.join(self.checkpoint_output, "info.txt"), 'r') as f:
            lines = f.readlines()

        module = self.current_step[0]
        step = self.current_step[1]
        lines[1] = f"{module},{step}\n"

        with open(os.path.join(self.checkpoint_output, "info.txt"), 'w') as f:
            f.writelines(lines)

    def module_step_execution_check(self, executing_module: int, executing_step: int) -> bool:
        return (self.checkpoint_input is None or
                self.current_step[0] < executing_module or
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

        if self.GRF_mode == "mix" and self.cpu_cores < 2 * prog_const.mix_short_seq_process_num:
            if self.GRF_mode == "smart":
                print("  \"native\" mode is selected due to insufficient number of available cpu cores.")
            else:
                print(f"   Number of available cpu cores insufficient "
                      f"(expect >= {2 * prog_const.mix_short_seq_process_num}"
                      f" but actually got {self.cpu_cores}), "
                      f"\"mix\" mode unavailable, redirect to \"native\" mode.")
            self.GRF_mode = "native"

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

        # Module 1, Step 1: Blast reference library in genome file
        if self.module_step_execution_check(1, 1):
            blast_reference.blast_genome_file(self)
            self.current_step = [1, 1]

        # Module 1, Step 2: Select 100% coverage entries from blast results
        if self.module_step_execution_check(1, 2):
            self.working_df_dict["base"] = process_homology.select_full_coverage(self)
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

        # Module 2, Step 1: Run TIRvish to find inverted repeats
        if self.module_step_execution_check(2, 1):
            print("Module 2, Step 1: Run TIRvish to find inverted repeats")
            self.working_df_dict["TIRvish"] = run_TIRvish.execute(self)
            self.current_step = [2, 1]
            self.save_checkpoint_file()

        # Module 2, Step 2: Process TIRvish results
        if self.module_step_execution_check(2, 2):
            print("Module 2, Step 2: Process TIRvish results")
            self.working_df_dict["TIRvish"] = process_de_novo_result.process_TIRvish_result(self)
            self.current_step = [2, 2]
            self.save_checkpoint_file()

        # Module 2, Step 3: Run GRF to find inverted repeats
        if self.module_step_execution_check(2, 3):
            print("Module 2, Step 3: Run GRF to find inverted repeats")
            self.working_df_dict["GRF"] = run_GRF.execute(self)
            self.current_step = [2, 3]
            self.save_checkpoint_file()

        # Module 2, Step 4: Process GRF results
        if self.module_step_execution_check(2, 4):
            print("Module 2, Step 4: Process GRF results")
            self.working_df_dict["GRF"] = process_de_novo_result.process_GRF_result(self)
            self.current_step = [2, 4]

        # Module 2, Step 5: Combine TIRvish and GRF results
        if self.module_step_execution_check(2, 5):
            print("Module 2, Step 5: Combine TIRvish and GRF results")
            process_de_novo_result.combine_de_novo_result(self)
            self.current_step = [2, 5]
            self.save_processed_de_novo_result_checkpoint_file()

        # Module 2, Step 6: Blast GRF and TIRvish result in reference library
        if self.module_step_execution_check(2, 6):
            # Checkpoint saving for this step is currently not available
            blast_reference.blast_de_novo_result(self)
            self.current_step = [2, 6]

        # Module 2, Step 7: Select 80% similar entries from blast results
        if self.module_step_execution_check(2, 7):
            self.working_df_dict["base"] = process_homology.select_eighty_similarity(self)
            self.current_step = [2, 7]
            self.save_checkpoint_file()

        # Module 2, Step 8: Get FASTA sequences from 80% similarity
        if self.module_step_execution_check(2, 8):
            print("Module 2, Step 8: Get FASTA sequences from 80% similarity")
            self.working_df_dict["base"] = get_fasta_sequence.execute(self)
            self.working_df_dict["m2_homo"] = self.working_df_dict["base"].copy()
            self.current_step = [2, 8]
            self.save_checkpoint_file()

        # Module 2, Step 6: Check TIR and TSD
        if self.module_step_execution_check(2, 9):
            print("Module 2, Step 9: Check TIR and TSD")
            self.working_df_dict["base"] = check_TIR_TSD.execute(self, module)
            self.current_step = [2, 9]
            self.save_checkpoint_file()

        # Module 2, Step 7: Save module result
        if self.module_step_execution_check(2, 10):
            print("Module 2, Step 10: Save module result")
            self.working_df_dict["m2"] = self.working_df_dict["base"]
            del self.working_df_dict["base"]
            self.current_step = [2, 10]
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

        # Module 3, Step 1: Prepare data
        if self.module_step_execution_check(3, 1):
            print("Module 3, Step 1: Prepare data")
            self.working_df_dict["base"] = prepare_data.execute(self, self.working_df_dict["m2_homo"])
            del self.working_df_dict["m2_homo"]
            self.current_step = [3, 1]
            self.save_checkpoint_file()

        # Module 3, Step 2: CNN prediction
        if self.module_step_execution_check(3, 2):
            print("Module 3, Step 2: CNN prediction")
            self.working_df_dict["base"] = CNN_predict.execute(self)
            self.current_step = [3, 2]
            self.save_checkpoint_file()

        # Module 3, Step 3: Get FASTA sequences from CNN prediction
        if self.module_step_execution_check(3, 3):
            print("Module 3, Step 3: Get FASTA sequences from CNN prediction")
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

        # Module 4, Step 1: Run TIRvish to find inverted repeats
        if self.module_step_execution_check(4, 1):
            print("Module 4, Step 1: Run TIRvish to find inverted repeats")
            self.working_df_dict["TIRvish"] = run_TIRvish.execute(self)
            self.current_step = [4, 1]
            self.save_checkpoint_file()

        # Module 4, Step 2: Process TIRvish results
        if self.module_step_execution_check(4, 2):
            print("Module 4, Step 2: Process TIRvish results")
            self.working_df_dict["TIRvish"] = process_de_novo_result.process_TIRvish_result(self)
            self.current_step = [4, 2]
            self.save_checkpoint_file()

        # Module 4, Step 3: Run GRF to find inverted repeats
        if self.module_step_execution_check(4, 3):
            print("Module 4, Step 3: Run GRF to find inverted repeats")
            self.working_df_dict["GRF"] = run_GRF.execute(self)
            self.current_step = [4, 3]
            self.save_checkpoint_file()

        # Module 4, Step 4: Process GRF results
        if self.module_step_execution_check(4, 4):
            print("Module 4, Step 4: Process GRF results")
            self.working_df_dict["GRF"] = process_de_novo_result.process_GRF_result(self)
            self.current_step = [4, 4]

        # Module 4, Step 5: Combine TIRvish and GRF results
        if self.module_step_execution_check(4, 5):
            print("Module 4, Step 5: Combine TIRvish and GRF results")
            process_de_novo_result.combine_de_novo_result(self)
            self.current_step = [4, 5]
            self.save_processed_de_novo_result_checkpoint_file()

        # Module 4, Step 6: Prepare data
        if self.module_step_execution_check(4, 6):
            print("Module 4, Step 6: Prepare data")
            self.working_df_dict["base"] = prepare_data.execute(self)
            self.current_step = [4, 6]
            self.save_checkpoint_file()

        # Module 4, Step 7: CNN prediction
        print("Module 4, Step 7: CNN prediction")
        if self.module_step_execution_check(4, 7):
            self.working_df_dict["base"] = CNN_predict.execute(self)
            self.current_step = [4, 7]
            self.save_checkpoint_file()

        # Module 4, Step 8: Get FASTA sequences from CNN prediction
        if self.module_step_execution_check(4, 8):
            print("Module 4, Step 8: Get FASTA sequences from CNN prediction")
            self.working_df_dict["base"] = get_fasta_sequence.execute(self)
            self.current_step = [4, 8]
            self.save_checkpoint_file()

        # Module 4, Step 9: Check TIR and TSD
        if self.module_step_execution_check(4, 9):
            print("Module 4, Step 9: Check TIR and TSD")
            self.working_df_dict["base"] = check_TIR_TSD.execute(self, module)
            self.current_step = [4, 9]
            self.save_checkpoint_file()

        # Module 4, Step 10: Save module result
        if self.module_step_execution_check(4, 10):
            self.working_df_dict["m4"] = self.working_df_dict["base"]
            del self.working_df_dict["base"]
            self.current_step = [4, 10]
            self.save_checkpoint_file()

        print("########################################################## Module 4 Finished "
              "######################################################")
        # return self.df_current.copy()
