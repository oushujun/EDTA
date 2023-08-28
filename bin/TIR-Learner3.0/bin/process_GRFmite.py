import os
import subprocess
import regex as re
import numpy as np
import pandas as pd
from Bio import SeqIO

spliter = "-+-"


def TA_repeats(s, percent=0.7):
    t = s.upper().count("T")
    a = s.upper().count("A")
    ta = t + a
    if (ta >= len(s) * percent):
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
    if (len(x["tsd"]) > 6 or x["tsd"] == "TAA" or x["tsd"] == "TTA"
        or x["tsd"] == "TA" or x["seq"][0:4] == "CACT") or x["seq"][0:4] == "GTGA":
        return True
    return np.nan


def execute(args):
    genome_name = args[1]

    # GRFmite_file_path = os.path.join("{genome_name}_GRFmite", "candidate.fasta")
    GRFmite_file_path = os.path.join(f"{genome_name}_filtered.fa_GRFmite", "candidate.fasta")

    df = pd.DataFrame({"id": [rec.id for rec in SeqIO.parse(GRFmite_file_path, "fasta")],
                       "seq": [str(rec.seq) for rec in SeqIO.parse(GRFmite_file_path, "fasta")],
                       "len": [len(rec) for rec in SeqIO.parse(GRFmite_file_path, "fasta")]})
    subprocess.Popen(["rm", "-rf", f"{genome_name}_filtered.fa_GRFmite"])
    df = df[df["len"] >= 50].copy()

    df["tir_len"] = df.swifter.progress_bar(True).apply(lambda x: find_digits_sum(x["id"].split(":")[-2]), axis=1)
    #df["tirLen"] = df["tirLen"].astype(int)
    df["tir"] = df.swifter.progress_bar(True).apply(lambda x: x["seq"][0:x["tir_len"]], axis=1)
    df["TA_repeats_seq_check"] = df.swifter.progress_bar(True).apply(lambda x:
                                                                    np.nan if TA_repeats(x["seq"]) else False, axis=1)
    df["check_N_per_seq_check"] = df.swifter.progress_bar(True).apply(lambda x:
                                                                    np.nan if check_N_per(x["seq"]) else False, axis=1)
    df["TA_repeats_tir_check"] = df.swifter.progress_bar(True).apply(lambda x:
                                                                    np.nan if TA_repeats(x["tir"]) else False, axis=1)
    df["check_N_tir_check"] = df.swifter.progress_bar(True).apply(lambda x:
                                                                 np.nan if check_N(x["tir"]) else False, axis=1)
    df = df.dropna(ignore_index=True).copy()

    df["tsd"] = df.swifter.progress_bar(True).apply(lambda x: x["id"].split(":")[-1], axis=1)
    df["tsd_check"] = df.swifter.progress_bar(True).apply(TSD_check, axis=1)
    df = df.dropna(ignore_index=True).loc[:, ["id", "seq"]].copy()

    df["id"] = df.swifter.progress_bar(True).apply(lambda x: ">" + x["id"], axis=1)
    # df.to_csv(os.path.join("../", f"{genome_name}{spliter}processedGRFmite.fa"), sep='\n', header=False, index=False)
    df.to_csv(f"{genome_name}{spliter}processedGRFmite.fa", sep='\n', header=False, index=False)
