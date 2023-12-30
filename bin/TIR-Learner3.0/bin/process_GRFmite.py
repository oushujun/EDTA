import os
import subprocess
import regex as re
import numpy as np
import pandas as pd
import swifter
from Bio import SeqIO

# from typing import TYPE_CHECKING
# if TYPE_CHECKING:
#     from main import TIRLearner

import prog_const


def TA_repeats(s, percent=0.7):
    t = s.upper().count("T")
    a = s.upper().count("A")
    ta = t + a
    if ta >= len(s) * percent:
        return True
    return False


def check_N(s):
    n = s.upper().count("N")
    if n > 0:
        return True
    else:
        return False


def check_N_per(s):
    n = s.upper().count("N")
    p = n / len(s)
    if p >= 0.20:
        return True
    else:
        return False


def find_digits_sum(string):
    pattern = '(\d+)'
    l = re.findall(pattern, string)
    return sum([int(i) for i in l])


def TSD_check(x):
    if ((len(x["tsd"]) > 6 or x["tsd"] == "TAA" or x["tsd"] == "TTA" or x["tsd"] == "TA" or x["seq"][0:4] == "CACT") or
            x["seq"][0:4] == "GTGA"):
        return True
    return np.nan


def execute(TIRLearner_instance):
    genome_name = TIRLearner_instance.genome_name
    flag_verbose = TIRLearner_instance.flag_verbose

    print("  Step 1/9: Retrieving GRFmite")
    # GRFmite_file_path = os.path.join("{genome_name}_GRFmite", "candidate.fasta")
    GRFmite_file_path = os.path.join(f"{genome_name}_filtered.fa_GRFmite", "candidate.fasta")

    df = pd.DataFrame({"id": [rec.id for rec in SeqIO.parse(GRFmite_file_path, "fasta")],
                       "seq": [str(rec.seq) for rec in SeqIO.parse(GRFmite_file_path, "fasta")],
                       "len": [len(rec) for rec in SeqIO.parse(GRFmite_file_path, "fasta")]})
    subprocess.Popen(["rm", "-rf", f"{genome_name}_filtered.fa_GRFmite"])
    df = df[df["len"] >= 50].copy()

    print("  Step 2/9: Getting TIR")
    df["tir_len"] = df.swifter.progress_bar(flag_verbose).apply(lambda x: find_digits_sum(x["id"].split(":")[-2]),
                                                                axis=1)
    # df["tirLen"] = df["tirLen"].astype(int)
    df["tir"] = df.swifter.progress_bar(flag_verbose).apply(lambda x: x["seq"][0:x["tir_len"]], axis=1)

    print("  Step 3/9: Checking TA repeats on sequence")
    df["TA_repeats_seq_check"] = df.swifter.progress_bar(flag_verbose).apply(lambda x:
                                                                             np.nan if TA_repeats(x["seq"]) else False,
                                                                             axis=1)
    print("  Step 4/9: Checking percentage of N on sequence")
    df["check_N_per_seq_check"] = df.swifter.progress_bar(flag_verbose).apply(lambda x:
                                                                              np.nan if check_N_per(
                                                                                  x["seq"]) else False, axis=1)
    print("  Step 5/9: Checking TA repeats on TIR")
    df["TA_repeats_tir_check"] = df.swifter.progress_bar(flag_verbose).apply(lambda x:
                                                                             np.nan if TA_repeats(x["tir"]) else False,
                                                                             axis=1)
    print("  Step 6/9: Checking N existance on TIR")
    df["check_N_tir_check"] = df.swifter.progress_bar(flag_verbose).apply(lambda x:
                                                                          np.nan if check_N(x["tir"]) else False,
                                                                          axis=1)
    df = df.dropna(ignore_index=True).copy()

    print("  Step 7/9: Getting TSD")
    df["tsd"] = df.swifter.progress_bar(flag_verbose).apply(lambda x: x["id"].split(":")[-1], axis=1)
    print("  Step 8/9: Checking TSD")
    df["tsd_check"] = df.swifter.progress_bar(flag_verbose).apply(TSD_check, axis=1)
    df = df.dropna(ignore_index=True).loc[:, ["id", "seq"]].copy()

    print("  Step 9/9: Saving processed GRFmite")
    df["id"] = df.swifter.progress_bar(flag_verbose).apply(lambda x: ">" + x["id"], axis=1)
    # df.to_csv(os.path.join("../", f"{genome_name}{spliter}processedGRFmite.fa"), sep='\n', header=False, index=False)
    df.to_csv(f"{genome_name}{prog_const.spliter}processedGRFmite.fa", sep='\n', header=False, index=False)
