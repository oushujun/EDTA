import sys
import os
import subprocess
import multiprocessing as mp
import pandas as pd

spliter = "-+-"
TIR_types = ("DTA", "DTC", "DTH", "DTM", "DTT")
blast_header = ("qacc", "sacc", "length", "pident", "gaps", "mismatch",
                "qstart", "qend", "sstart", "send", "evalue", "qcovhsp")
blast_type = {"length": int, "gaps": int, "mismatch": int,
                "qstart": int, "qend": int, "sstart": int, "send": int}


# def processHomology(genome_name, species):
#     homo = pd.DataFrame()
#     for i in TEs:
#         blast = f"{genome_name}{spliter}blast{spliter}{species}_{i}_RefLib"
#         if os.path.exists(blast) and os.path.getsize(blast) != 0:
#             df = pd.read_csv(blast, header=None, names=blast_header, dtype=blast_type, sep="\t")
#             df["sacc"] = df["sacc"].astype(str)
#             df = df.loc[(df["qcovhsp"] == 100) & (df["pident"] >= 80)]
#             df = df.sort_values(["sacc", "sstart", "send", "qcovhsp", "pident"],
#                                       ascending=[True, True, True, True, True])
#             df = df.drop_duplicates(["sacc", "sstart", "send"], keep="last")
#             df.insert(0, "TE", i)
#             homo = pd.concat([homo, df], ignore_index=True)
#     return homo


def process_homology(genome_name, species, TIR_type):
    blast = f"{genome_name}{spliter}blast{spliter}{species}_{TIR_type}_RefLib"
    df = None
    if os.path.exists(blast) and os.path.getsize(blast) != 0:
        # df = pd.read_csv(blast, sep='\t', header=None, names=blast_header, dtype=blast_type, engine="pyarrow")
        df = pd.read_csv(blast, sep='\t', header=None, names=blast_header, dtype=blast_type, engine='c',
                         memory_map=True)
        # df["sacc"] = df["sacc"].astype(str)
        df = df.loc[(df["qcovhsp"] == 100) & (df["pident"] >= 80)].reset_index(drop=True)
        df = df.sort_values(["sacc", "sstart", "send", "qcovhsp", "pident"],
                            ascending=[True, True, True, True, True], ignore_index=True)
        df = df.drop_duplicates(["sacc", "sstart", "send"], keep="last", ignore_index=True)
        df.insert(0, "TIR_type", TIR_type)
    return df


# def verifyResult(df_in, species):
#     if df_in.shape[0] == 0:
#         print(f"""
#         ERROR: No sequence is found similar to the TIR database of {species}!
#         You may have specified the wrong species. Please double check or set species=others and rerun TIR-Learner
#         """)
#         sys.exit(-1)

# def verifyResult(df_list, species):
#     if df_list == [None] * 5:
#         print(f"""
#         ERROR: No sequence is found similar to the TIR database of {species}!
#         You may have specified the wrong species. Please double check or set species=others and rerun TIR-Learner
#         """)
#         sys.exit(-1)
#
#
# def processResult(df_in):
#     df = df_in.iloc[:, [0, 1, 2, 9, 10]].copy()
#     df = df.set_axis(["TIR_type", "id", "seqid", "sstart", "send"], axis=1)
#     return df


def process_result(df_list, species):
    try:
        df = pd.concat(df_list, ignore_index=True).iloc[:, [0, 1, 2, 9, 10]].copy()
    except ValueError:
        print(f"""
        ERROR: No sequence is found similar to the TIR database of {species}!
        You may have specified the wrong species. Please double check or set species=others and rerun TIR-Learner. 
        """)
        sys.exit(-1)
    df = df.set_axis(["TIR_type", "id", "seqid", "sstart", "send"], axis=1)
    return df


def execute(args):
    print("Module 1, Step 2: Select 100% coverage entries from Blast results")
    genome_name = args[1]
    t = args[3]
    species = args[4]

    mp_args_list = [(genome_name, species, TIR_type) for TIR_type in TIR_types]
    with mp.Pool(int(t)) as pool:
        df_list = pool.starmap(process_homology, mp_args_list)
    # subprocess.Popen(["rm", "-f", f"*{spliter}blast{spliter}*"])  # remove blast files
    subprocess.Popen(["find", ".", "-name", f"*{spliter}blast{spliter}*", "-delete"])  # remove blast files
    # subprocess.Popen(f"rm -f *{spliter}blast{spliter}*", shell=True)
    return process_result(df_list, species)
