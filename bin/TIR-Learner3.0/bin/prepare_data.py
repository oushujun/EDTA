import subprocess
import pandas as pd
from Bio import SeqIO

# from typing import TYPE_CHECKING
# if TYPE_CHECKING:
#     from main import TIRLearner


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
    processedGRFmite_file = TIRLearner_instance.processedGRFmite_file

    data = [{"id": rec.id, "seq": rec.seq} for rec in SeqIO.parse(processedGRFmite_file, "fasta")]
    subprocess.Popen(["unlink", processedGRFmite_file])

    # Get non_homo data
    if df_homo is not None:
        data = [rec for rec in data if rec["id"] not in df_homo["id"].unique()]

    df = pd.DataFrame(data, columns=["id", "seq"], dtype=str)
    return df
