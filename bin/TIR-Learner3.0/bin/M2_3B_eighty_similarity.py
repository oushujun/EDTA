import os
import subprocess
import multiprocessing as mp
import pandas as pd

from M1_2_full_coverage import process_result

spliter = "-+-"
TIR_types = ("DTA", "DTC", "DTH", "DTM", "DTT")

blast_header = ("qseqid", "sseqid", "length", "pident", "gaps", "mismatch",
                "qstart", "qend", "sstart", "send", "evalue", "qcovhsp")
blast_type = {"length": int, "gaps": int, "mismatch": int,
                "qstart": int, "qend": int, "sstart": int, "send": int}


# def processHomology(file, species):
#     homo = pd.DataFrame()
#     for i in TEs:
#         blast = f"{file}{spliter}blast{spliter}{species}_{i}_RefLib"
#         if os.path.exists(blast) and os.path.getsize(blast) != 0:
#             df = pd.read_csv(blast, header=None, names=blast_header, dtype=blast_type, sep="\t")
#             df = df.loc[(df["qcovhsp"] >= 80) & (df["pident"] >= 80)]
#
#             df["sseqid"] = df.swifter.progress_bar(True).apply(lambda x: x["qseqid"].split(":")[0], axis=1)
#             df["sstart"] = df.swifter.progress_bar(True).apply(lambda x: int(x["qseqid"].split(":")[1]), axis=1)
#             df["send"] = df.swifter.progress_bar(True).apply(lambda x: int(x["qseqid"].split(":")[2]), axis=1)
#
#             df = df.sort_values(["sseqid", "sstart", "send", "qcovhsp", "pident"],
#                                 ascending=[True, True, True, True, True])
#             df = df.drop_duplicates(["sseqid", "sstart", "send"], keep="last")
#             df.insert(0, "TE", i)
#             homo = pd.concat([homo, df], ignore_index=True)
#     return homo


def process_homology(file_name, species, TIR_type):
    blast = f"{file_name}{spliter}blast{spliter}{species}_{TIR_type}_RefLib"
    df = None
    if os.path.exists(blast) and os.path.getsize(blast) != 0:
        # df = pd.read_csv(blast, sep='\t', header=None, names=blast_header, dtype=blast_type, engine="pyarrow")
        df = pd.read_csv(blast, sep='\t', header=None, names=blast_header, dtype=blast_type, engine='c',
                         memory_map=True)
        df = df.loc[(df["qcovhsp"] >= 80) & (df["pident"] >= 80)].reset_index(drop=True)

        df["sseqid"] = df.swifter.progress_bar(True).apply(lambda x: x["qseqid"].split(":")[0], axis=1)
        df["sstart"] = df.swifter.progress_bar(True).apply(lambda x: int(x["qseqid"].split(":")[1]), axis=1)
        df["send"] = df.swifter.progress_bar(True).apply(lambda x: int(x["qseqid"].split(":")[2]), axis=1)

        df = df.sort_values(["sseqid", "sstart", "send", "qcovhsp", "pident"],
                            ascending=[True, True, True, True, True], ignore_index=True)
        df = df.drop_duplicates(["sseqid", "sstart", "send"], keep="last", ignore_index=True)
        df.insert(0, "TIR_type", TIR_type)
    return df


def execute(args):
    print("Module 2, Step 3B: Select 80% similar entries from Blast results")
    genome_name = args[1]
    processedGRFmite_file_name = f"{genome_name}{spliter}processedGRFmite.fa"
    t = args[3]
    species = args[4]

    mp_args_list = [(processedGRFmite_file_name, species, TIR_type) for TIR_type in TIR_types]
    with mp.Pool(int(t)) as pool:
        df_list = pool.starmap(process_homology, mp_args_list)
    subprocess.Popen(["find", ".", "-name", f"*{spliter}blast{spliter}*", "-delete"])  # remove blast files
    return process_result(df_list, species)
