import subprocess
import pandas as pd
import swifter
from Bio import SeqIO

# from typing import TYPE_CHECKING
# if TYPE_CHECKING:
#     from main import TIRLearner


def get_start_end(genome_file, df_in, flag_verbose, length=200):
    df = df_in.copy()
    df["start"] = df.swifter.progress_bar(flag_verbose).apply(lambda x: min(x["sstart"], x["send"]), axis=1)
    df["end"] = df.swifter.progress_bar(flag_verbose).apply(lambda x: max(x["sstart"], x["send"]), axis=1)
    df["sstart"] = df["start"]
    df["send"] = df["end"]

    df["start"] = df.swifter.progress_bar(flag_verbose).apply(lambda x: max(x["sstart"] - length, 1) - 1, axis=1)
    # start = start - length if start > length else 1, start--
    df["end"] = df.swifter.progress_bar(flag_verbose).apply(lambda x: x["send"] + length, axis=1)

    # Ensure "end" not exceeding seq's length
    fasta_len_dict = {rec.id: len(rec.seq) for rec in SeqIO.parse(genome_file, "fasta")}
    df["end"] = df.swifter.progress_bar(flag_verbose).apply(lambda x: min(x["end"], fasta_len_dict[x["seqid"]]), axis=1)

    # df["p_start"] = df.apply(lambda x: max(min(x["p1"], x["p2"]) - length, 1) - 1, axis=1)
    # # p_start = min(p1, p2)
    # # p_start = p_start - 200 if p_start > 200 else 1
    # df["p_end"] = df.apply(lambda x: max(x["p1"], x["p2"]) + length, axis=1)

    df = df.drop_duplicates(["start", "end", "seqid", "TIR_type"], keep="first", ignore_index=True)
    return df


# def getFastaPieces(genome_file, df_in):
#     df = df_in.reset_index(drop=True)
#     # df_bed = df.iloc[:, [2, 5, 6]]
#     df_bed = df.loc[:, ["seqid", "start", "end"]]
#     bed_string = df_bed.to_string(header=False, index=False)
#     bed = BedTool(bed_string, from_string=True)
#     fasta = BedTool(genome_file)
#     bed = bed.sequence(fi=fasta)
#     bed_seq_file = bed.seqfn
#     bed_series = pd.read_csv(bed_seq_file, header=None).squeeze()
#     df["seq"] = bed_series[1::2].reset_index(drop=True)
#     return df


def get_fasta_pieces_bedtools(genome_file, df_in: pd.DataFrame):
    df = df_in.reset_index(drop=True)
    # df_bed = df.iloc[:, [2, 5, 6]]
    df.loc[:, ["seqid", "start", "end"]].to_csv("bed.txt", sep='\t', header=False, index=False)
    # subprocess.Popen(f"bedtools getfasta -fo seq_from_bed.txt -fi {genome_file} -bed bed.txt", ).wait()
    subprocess.Popen(["bedtools", "getfasta", "-fo", "seq_from_bed.txt", "-fi", genome_file, "-bed", "bed.txt"]).wait()
    # bed_series = pd.read_csv("seq_from_bed.txt", header=None, engine="pyarrow").squeeze()
    bed_series = pd.read_csv("seq_from_bed.txt", header=None, engine='c', memory_map=True).squeeze()
    subprocess.Popen(["find", ".", "-name", "*bed.txt", "-delete"])
    df["seq"] = bed_series[1::2].reset_index(drop=True)
    return df


def execute(TIRLearner_instance) -> pd.DataFrame:
    df_in = TIRLearner_instance.working_df_dict["base"]
    genome_file = TIRLearner_instance.genome_file
    flag_verbose = TIRLearner_instance.flag_verbose

    df = get_start_end(genome_file, df_in, flag_verbose)
    return get_fasta_pieces_bedtools(genome_file, df)
