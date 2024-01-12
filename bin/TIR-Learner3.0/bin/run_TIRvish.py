import os
import subprocess
import pandas as pd
import swifter  # ATTENTION: DO NOT REMOVE "swifter" EVEN IF IDE SHOWS IT IS NOT USED!
from get_fasta_sequence import get_fasta_pieces_SeqIO

import prog_const


def run_TIRvish(genome_file, genome_name, TIR_length, gt_path):
    gt_bin_path = os.path.join(gt_path, "gt")
    gt_index_file_name = genome_name + prog_const.spliter + "gt_index"
    subprocess.Popen(
        [gt_bin_path, "suffixerator", "-db", genome_file, "-indexname", gt_index_file_name,
         "-tis", "-suf", "-lcp", "-des", "-ssp", "-sds", "-dna", "-mirrored"]).wait()

    TIRvish_result_gff3_file_name = f"{genome_name}{prog_const.spliter}TIRvish.gff3"
    gt_tirvish = (f"\"{gt_bin_path}\" tirvish -index {gt_index_file_name} -seed 20 -mintirlen 10 -maxtirlen 1000 "
                  f"-mintirdist 10 -maxtirdist {str(TIR_length)} -similar 80 -mintsd 2 -maxtsd 11 "
                  f"-vic 13 -seqids \"yes\" > {TIRvish_result_gff3_file_name}")
    subprocess.Popen(gt_tirvish, shell=True).wait()
    subprocess.Popen(["find", ".", "-name", f"{gt_index_file_name}*", "-delete"])
    return TIRvish_result_gff3_file_name


def get_TIRvish_result_df(TIRvish_result_gff3_file_name):
    df_data_dict = {"seqid": [], "start": [], "end": [], "TIR1_start": [], "TIR1_end": [], "TIR2_start": [],
                    "TIR2_end": [], "id": []}
    df_type = {"start": int, "end": int, "TIR1_start": int, "TIR1_end": int, "TIR2_start": int, "TIR2_end": int}

    if os.path.exists(TIRvish_result_gff3_file_name) and os.path.getsize(TIRvish_result_gff3_file_name) != 0:
        with open(TIRvish_result_gff3_file_name, 'r') as f:
            while (line := f.readline()) != "":
                while line.startswith('#'):
                    line = f.readline()
                TSD1 = list(map(int, f.readline().split('\t')[3:5]))
                record = f.readline().split('\t')
                TIR1 = list(map(int, f.readline().split('\t')[3:5]))
                TIR2 = list(map(int, f.readline().split('\t')[3:5]))
                TSD2 = list(map(int, f.readline().split('\t')[3:5]))
                f.readline()  # Jump the line "###"

                seqid = record[0]
                start = int(record[3])
                end = int(record[4])
                # info = record[8]
                # df_data_dict["info"].append(info)
                df_data_dict["seqid"].append(seqid)
                df_data_dict["start"].append(start - 1)
                df_data_dict["end"].append(end)
                df_data_dict["TIR1_start"].append(TIR1[0] - 1)
                df_data_dict["TIR1_end"].append(TIR1[1])
                df_data_dict["TIR2_start"].append(TIR2[0] - 1)
                df_data_dict["TIR2_end"].append(TIR2[1])
                df_data_dict["id"].append(
                    f">{seqid}:{start}:{end}:"
                    f"tsd1_{TSD1[0]}_{TSD1[1]}_tsd2_{TSD2[0]}_{TSD2[1]}_"
                    f"tir1_{TIR1[0]}_{TIR1[1]}_tir2_{TIR2[0]}_{TIR2[1]}")

    subprocess.Popen(["unlink", TIRvish_result_gff3_file_name])
    return pd.DataFrame(df_data_dict).astype(df_type)


def execute(TIRLearner_instance):
    genome_file = TIRLearner_instance.genome_file_path
    genome_name = TIRLearner_instance.genome_name
    cpu_cores = TIRLearner_instance.cpu_cores
    TIR_length = TIRLearner_instance.TIR_length
    flag_verbose = TIRLearner_instance.flag_verbose
    gt_path = TIRLearner_instance.gt_path

    print("  Step 1/2: Executing TIRvish")
    TIRvish_result_gff3_file_name = run_TIRvish(genome_file, genome_name, TIR_length, gt_path)
    print("  Step 2/2: Getting TIRvish result")
    df = get_TIRvish_result_df(TIRvish_result_gff3_file_name)
    return get_fasta_pieces_SeqIO(genome_file, df, cpu_cores, flag_verbose)
