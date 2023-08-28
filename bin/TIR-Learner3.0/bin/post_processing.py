import os
import multiprocessing as mp
import numpy as np
import pandas as pd

from get_fasta_sequence import get_fasta_pieces_bedtools

spliter = "-+-"


def TA_repeats(TIR_pair, percent=0.7):
    s = TIR_pair[0] + TIR_pair[1]
    t = s.upper().count("T")
    a = s.upper().count("A")
    ta = t + a
    if (ta >= len(s) * percent):
        return True
    return False


def combine_all(df_list):
    if len(df_list) > 1:
        df = pd.concat(df_list, ignore_index=True)
    else:
        df = df_list[0].copy()
    df = df.sort_values(["seqid", "sstart", "send", "source", "type"], ignore_index=True)
    df = df.drop_duplicates(["seqid", "sstart", "send"], keep="first", ignore_index=True)
    df = df.drop_duplicates(["seqid", "sstart"], keep="first", ignore_index=True)
    df = df.drop_duplicates(["seqid", "send"], keep="first", ignore_index=True)
    df["TArepeats_TIR_check"] = df.swifter.progress_bar(True).apply(lambda x: np.nan if TA_repeats(x["TIR"]) else False,
                                                                    axis=1)
    df = df.dropna(ignore_index=True).drop(columns="TArepeats_TIR_check")
    return df


def format_dataframe_in_gff3_format(df_in):
    df = df_in.copy()
    df["attributes"] = df.swifter.progress_bar(True).apply(lambda x: f"TIR:{x['TIR'][0]}_{x['TIR'][1]}_{x['p_TIR']}_" +
                                                                     f"TSD:{x['TSD'][0]}_{x['TSD'][1]}_{x['p_TSD']}" +
                                                                     f"{spliter}{x['len']}", axis=1)
    df = df.loc[:, ["seqid", "source", "type", "sstart", "send", "attributes"]]
    df.insert(5, "phase", ".")
    df.insert(5, "strand", ".")
    df.insert(5, "score", ".")
    return df


# =============================== Remove the Shorter One Among Two Overlapped Sequences ================================

def check_element_overlap(x1, y1, x2, y2):
    # For only elements
    # Precondition: x1 < x2, y1 != y2
    # True: Overlap, x2 <= y1 < y2
    # x1------A------y1
    #        x2------B------y2
    # False: Nesting, y2 < y1
    # x1--------A---------y1
    #        x2----B----y2
    # False: No correlation, y1 < x2
    # x1----A----y1
    #                x2------B------y2
    if y1 >= x2 and y1 < y2:
        return True
    return False


def check_element_TIR_overlap(x1, y1, x2, y2, m1, n1, m2, n2):
    # For elements and TIRs
    # Precondition: x1 < x2, y1 < x2 or y1 > y2
    # True: Overlap - A right element with B left TIR
    # m1===x1------A------y1===n1
    #                  m2========x2------B------y2========n2
    # True: Overlap - B left element with A right TIR
    # m1========x1------A------y1========n1
    #                            m2===x2------B------y2===n2
    # True: Overlap - A right element with B right TIR
    # m1===x1----------------A----------------y1===n1
    #             m2========x2------B------y2========n2
    # True: Overlap - A left element with B left TIR
    #   m1===x1----------------A----------------y1===n1
    # m2========x2------B------y2========n2
    if y1 < x2:
        if y1 > m2 or x2 < n1:
            return True
    elif y1 > y2:
        if y1 < n2 or x1 > m2:
            return True
    return False


def remove_overlap(df_in):
    df = df_in.sort_values(by=["sstart", "send"], ignore_index=True)
    df["TIR_len"] = df.swifter.progress_bar(True).apply(lambda x: len(x["TIR"][0]), axis=1)
    df["tstart"] = df.loc[:, "sstart"] - df.loc[:, "TIR_len"]
    df["tend"] = df.loc[:, "send"] + df.loc[:, "TIR_len"]
    idx_sstart = df.columns.get_loc("sstart")
    idx_send = df.columns.get_loc("send")
    idx_tstart = df.columns.get_loc("tstart")
    idx_tend = df.columns.get_loc("tend")
    ptr = 0
    while ptr+1 < df.shape[0]:
        if (check_element_overlap(df.iloc[ptr, idx_sstart], df.iloc[ptr, idx_send],
                                  df.iloc[ptr+1, idx_sstart], df.iloc[ptr+1, idx_send]) or
                check_element_TIR_overlap(df.iloc[ptr, idx_sstart], df.iloc[ptr, idx_send],
                                          df.iloc[ptr+1, idx_sstart], df.iloc[ptr+1, idx_send],
                                          df.iloc[ptr, idx_tstart], df.iloc[ptr, idx_tend],
                                          df.iloc[ptr+1, idx_tstart], df.iloc[ptr+1, idx_tend])):
            df = df.drop(labels=df.iloc[[ptr, ptr+1], df.columns.get_loc("len")].idxmin()).reset_index(drop=True)
        ptr += 1
    df = df.drop(columns=["TIR_len", "tstart", "tend"])
    return df

# ======================================================================================================================


def get_final_fasta_file(df_in, genome_file, genome_name, file_path):
    df = df_in.copy()
    df["name"] = df.swifter.progress_bar(True).apply(lambda x: (f">{genome_name}_{x['seqid']}_{x['sstart']}_{x['send']}"
                                                                f"_{x['type']}_{x['attributes']}"), axis=1)
    df.rename(columns={"sstart": "start", "send": "end"}, inplace=True)
    df = get_fasta_pieces_bedtools(genome_file, df)
    df = df.loc[:, ["name", "seq"]]
    df.to_csv(file_path, index=False, header=False, sep="\n")


def execute(args, df_list, output_dir):
    genome_file = args[0]
    genome_name = args[1]
    t = args[3]

    print("############################################################ Post Processing  "
          "#########################################################")
    output_folder_path = os.path.join(output_dir, "TIR-Learner-Result")
    os.makedirs(output_folder_path, exist_ok=True)

    df_combined = combine_all(df_list)
    df_gff3 = format_dataframe_in_gff3_format(df_combined)
    df_gff3.to_csv(os.path.join(output_folder_path, f"{genome_name}_FinalAnn.gff3"), index=False, header=False,
                   sep="\t")
    get_final_fasta_file(df_gff3, genome_file, genome_name, os.path.join(output_folder_path,
                                                                         f"{genome_name}_FinalAnn.fa"))
    del df_gff3

    df_combined_groupby_seqid = df_combined.groupby("seqid")
    df_combined_seqid_list = [df_combined_groupby_seqid.get_group(df) for df in df_combined_groupby_seqid.groups]
    del df_combined, df_combined_groupby_seqid
    with mp.Pool(int(t)) as pool:
        df_filtered_list = pool.map(remove_overlap, df_combined_seqid_list)
    df_filtered = pd.concat(df_filtered_list, ignore_index=True)
    del df_combined_seqid_list, df_filtered_list

    df_gff3_filtered = format_dataframe_in_gff3_format(df_filtered)
    df_gff3_filtered.to_csv(os.path.join(output_folder_path, f"{genome_name}_FinalAnn_filter.gff3"), index=False,
                            header=False, sep="\t")
    get_final_fasta_file(df_gff3_filtered, genome_file, genome_name,
                         os.path.join(output_folder_path, f"{genome_name}_FinalAnn_filter.fa"))
    del df_gff3_filtered
