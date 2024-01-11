import os
import subprocess
import multiprocessing as mp
import pandas as pd

import prog_const

blast_header_full_coverage = ("qacc", "sacc", "length", "pident", "gaps", "mismatch",
                              "qstart", "qend", "sstart", "send", "evalue", "qcovhsp")

blast_header_eighty_similarity = ("qseqid", "sseqid", "length", "pident", "gaps", "mismatch",
                                  "qstart", "qend", "sstart", "send", "evalue", "qcovhsp")

blast_type = {"length": int, "gaps": int, "mismatch": int,
              "qstart": int, "qend": int, "sstart": int, "send": int}


def process_homology_full_coverage(genome_name, species, TIR_type):
    blast = f"{genome_name}{prog_const.spliter}blast{prog_const.spliter}{species}_{TIR_type}_RefLib"
    df = None
    if os.path.exists(blast) and os.path.getsize(blast) != 0:
        # df = pd.read_csv(blast, sep='\t', header=None, names=blast_header_full_coverage, dtype=blast_type, engine="pyarrow")
        df = pd.read_csv(blast, sep='\t', header=None, names=blast_header_full_coverage, dtype=blast_type, engine='c',
                         memory_map=True)
        # df["sacc"] = df["sacc"].astype(str)
        df = df.loc[(df["qcovhsp"] == 100) & (df["pident"] >= 80)].reset_index(drop=True)
        df = df.sort_values(["sacc", "sstart", "send", "qcovhsp", "pident"],
                            ascending=[True, True, True, True, True], ignore_index=True)
        df = df.drop_duplicates(["sacc", "sstart", "send"], keep="last", ignore_index=True)
        df.insert(0, "TIR_type", TIR_type)
    return df


def process_homology_eighty_similarity(file_name, species, TIR_type):
    blast = f"{file_name}{prog_const.spliter}blast{prog_const.spliter}{species}_{TIR_type}_RefLib"
    df = None
    if os.path.exists(blast) and os.path.getsize(blast) != 0:
        # df = pd.read_csv(blast, sep='\t', header=None, names=blast_header_eighty_similarity, dtype=blast_type, engine="pyarrow")
        df = pd.read_csv(blast, sep='\t', header=None, names=blast_header_eighty_similarity, dtype=blast_type,
                         engine='c', memory_map=True)
        df = df.loc[(df["qcovhsp"] >= 80) & (df["pident"] >= 80)].reset_index(drop=True)

        df["sseqid"] = df.swifter.progress_bar(True).apply(lambda x: x["qseqid"].split(":")[0], axis=1)
        df["sstart"] = df.swifter.progress_bar(True).apply(lambda x: int(x["qseqid"].split(":")[1]), axis=1)
        df["send"] = df.swifter.progress_bar(True).apply(lambda x: int(x["qseqid"].split(":")[2]), axis=1)

        df = df.sort_values(["sseqid", "sstart", "send", "qcovhsp", "pident"],
                            ascending=[True, True, True, True, True], ignore_index=True)
        df = df.drop_duplicates(["sseqid", "sstart", "send"], keep="last", ignore_index=True)
        df.insert(0, "TIR_type", TIR_type)
    return df


def process_result(df_list, species):
    try:
        df = pd.concat(df_list, ignore_index=True).iloc[:, [0, 1, 2, 9, 10]].copy()
    except ValueError:
        # print(f"""
        # ERROR: No sequence is found similar to the TIR database of {species}!
        # You may have specified the wrong species. Please double-check or set species=others and rerun TIR-Learner.
        # """)
        # sys.exit(-1)
        raise SystemExit(f"""
        ERROR: No sequence is found similar to the TIR database of {species}!
        You may have specified the wrong species. Please double-check or set species=others and rerun TIR-Learner. 
        """)
    df = df.set_axis(["TIR_type", "id", "seqid", "sstart", "send"], axis=1)
    return df


def select_full_coverage(TIRLearner_instance) -> pd.DataFrame:
    print("Module 1, Step 2: Select 100% coverage entries from blast results")
    mp_args_list = [(TIRLearner_instance.genome_name, TIRLearner_instance.species, TIR_type)
                    for TIR_type in prog_const.TIR_types]
    with mp.Pool(int(TIRLearner_instance.cpu_cores)) as pool:
        df_list = pool.starmap(process_homology_full_coverage, mp_args_list)
    # subprocess.Popen(["rm", "-f", f"*{spliter}blast{spliter}*"])  # remove blast files
    subprocess.Popen(["find", ".", "-name", f"*{prog_const.spliter}blast{prog_const.spliter}*", "-delete"])
    # subprocess.Popen(f"rm -f *{spliter}blast{spliter}*", shell=True)
    return process_result(df_list, TIRLearner_instance.species)


def select_eighty_similarity(TIRLearner_instance) -> pd.DataFrame:
    print("Module 2, Step 7: Select 80% similar entries from blast results")
    mp_args_list = [(TIRLearner_instance.processed_de_novo_result_file_name, TIRLearner_instance.species, TIR_type)
                    for TIR_type in prog_const.TIR_types]
    with mp.Pool(int(TIRLearner_instance.cpu_cores)) as pool:
        df_list = pool.starmap(process_homology_eighty_similarity, mp_args_list)
    subprocess.Popen(["find", ".", "-name", f"*{prog_const.spliter}blast{prog_const.spliter}*", "-delete"])
    return process_result(df_list, TIRLearner_instance.species)
