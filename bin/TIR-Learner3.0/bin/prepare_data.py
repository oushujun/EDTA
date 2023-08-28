import subprocess
import pandas as pd
from Bio import SeqIO

spliter = "-+-"

# def getNonHomo(file, df_homo):
#     non_homo = []
#     for rec in SeqIO.parse(file, "fasta"):
#         if rec.id not in df_homo["genome"].unique():
#             non_homo.append({"genome": rec.id, "seq": rec.seq})
#     df_non_homo = pd.DataFrame(non_homo, columns=["genome", "seq"], dtype=str)
#     return df_non_homo
#
#
# def getAll(file):
#     all = []
#     for rec in SeqIO.parse(file, "fasta"):
#             all.append({"genome": rec.id, "seq": rec.seq})
#     df_all = pd.DataFrame(all, columns=["genome", "seq"], dtype=str)
#     return df_all


def execute(args, df_homo=None):
    genome_name = args[1]
    # processedGRFmite_file_dir = os.path.join("../", f"{genome_name}{spliter}processedGRFmite.fa")
    processedGRFmite_file_name = f"{genome_name}{spliter}processedGRFmite.fa"

    data = [{"id": rec.id, "seq": rec.seq} for rec in SeqIO.parse(processedGRFmite_file_name, "fasta")]
    subprocess.Popen(["unlink", processedGRFmite_file_name])

    # Get non_homo data
    if df_homo is not None:
        data = [rec for rec in data if rec["id"] not in df_homo["id"].unique()]

    df = pd.DataFrame(data, columns=["id", "seq"], dtype=str)
    return df
