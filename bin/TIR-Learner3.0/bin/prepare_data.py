import subprocess
import pandas as pd
from Bio import SeqIO


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


def execute(TIRLearner_instance, df_homo: pd.DataFrame = None) -> pd.DataFrame:
    data = [{"id": rec.id, "seq": rec.seq}
            for rec in SeqIO.parse(TIRLearner_instance.processed_de_novo_result_file_name, "fasta")]
    subprocess.Popen(["unlink", TIRLearner_instance.processed_de_novo_result_file_name])

    # Get non_homo data
    if df_homo is not None:
        data = [rec for rec in data if rec["id"] not in df_homo["id"].unique()]

    df = pd.DataFrame(data, columns=["id", "seq"], dtype=str)
    return df
