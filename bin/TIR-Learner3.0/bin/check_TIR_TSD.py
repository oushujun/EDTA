import re
import numpy as np
from Bio.Seq import Seq

spliter = "-+-"

#
TSD = {}
TSD["DTA"] = [8]
TSD["DTC"] = [3, 2]
TSD["DTH"] = [3]
TSD["DTM"] = [10, 9, 8, 7]
TSD["DTT"] = [2]

# DTA:8
# DTC:2/3
# DTH:3(twa)
# DTM:7-10
# DTT:2(TA)

def compare(tir1, tir2):
    d = 0
    for i in range(0, len(tir1)):
        if (tir1[i] != tir2[i]):
            d += 1
    return d


def sliding_window(seq1, seq2, tsdlength):
    set1 = []
    set2 = []
    for i in range(0, len(seq1) - tsdlength + 1):
        set1.append(seq1[i:i + tsdlength])
        set2.append(seq2[i:i + tsdlength])
    return set1, set2


def conserved(fam, s1):
    noMotif = ["DTX", "NonTIR"]
    if (fam in noMotif):
        return True
    else:
        motif1 = " "
        pattern = " "
        if fam == "DTA":
            # YARNG
            motif1 = s1.upper()[0:5]
            pattern = "[CT]A[AG][ATGC]G"
        if fam == "DTC":
            # CMCWR
            motif1 = s1.upper()[0:5]
            pattern = "CACT[AG]"
        if fam == "DTH":
            motif1 = s1.upper()[0:4]
            pattern = "G[GA][GC]C"
        if fam == "DTM":
            motif1 = s1.upper()[0:1]
            pattern = "[GC]"
        if fam == "DTT":
            motif1 = s1.upper()[0:10]
            pattern = "CT[ATCG][ATCG]CTC[ATCG][ATCG]T"
        if fam == "DTE":
            # GGNRM
            motif1 = s1.upper()[0:5]
            pattern = "GG[ATCG][AG][AC]"
        if fam == "DTR":
            # CACWATG
            motif1 = s1.upper()[0:7]
            pattern = "CAC[AT]ATG"
        if fam == "DTP":
            # CANRG
            motif1 = s1.upper()[0:5]
            pattern = "CA[ATGC][AG]G"
        motif2 = str(Seq(motif1).reverse_complement())
        z1 = bool(re.match(pattern, motif1))
        z2 = bool(re.match(pattern, motif2))
        if (z1 == True or z2 == True):
            return True
        else:
            return False


def conserved_DTH(set1, tsd_dffset, l):
    for i in tsd_dffset:
        if (tsd_dffset[i] < l * 0.2):
            s1 = set1[int(i.split(":")[0])]
            if (s1.lower() == "tta" or s1.lower == "taa"):
                return True
    return False


def conserved_DTT(set1, tsd_dffset, l):
    for i in tsd_dffset:
        if (tsd_dffset[i] < l * 0.2):
            s1 = set1[int(i.split(":")[0])]
            if (s1[0:2].lower() == "ta"):
                return True
    return False


def get_difference(set1, set2):
    tsd_diff = {}
    for i in range(0, len(set1)):
        for j in range(0, len(set2)):
            name = str(i) + ":" + str(j)
            diff = compare(set1[i], set2[j])
            tsd_diff[name] = diff
    return tsd_diff


def is_TSD(tsd_dffset, l):
    for i in tsd_dffset:
        if tsd_dffset[i] < l * 0.2:
            return True
    return False


def check_TIR(x):
    family = x[0]
    s = x["seq"][200:-200].upper()
    len_s = len(s)
    minL = 10

    if (len_s <= 200):
        l_List = list(range(minL, int(len_s / 2)))
    else:
        l_List = list(range(minL, 100))

    for l in l_List:
        s1 = s[0:l]
        s2_ = s[-l:]
        s2 = Seq(s2_).reverse_complement()
        d = compare(s1, s2)
        if d < l * 0.2 and conserved(family, s1):
            return l
    return np.nan


def check_TSD(x):
    family = x[0]
    s = x["seq"]
    l = TSD[family]

    for i in l:
        s1 = s[200 - i:200]
        last20 = s[-200:]
        s2 = last20[0:i]
        set1, set2 = sliding_window(s1, s2, i)
        dff = get_difference(set1, set2)
        if family == "DTH" and conserved_DTH(set1, dff, i):
            return i
        elif family == "DTT" and conserved_DTT(set1, dff, i):
            return i
        elif family != "DTH" and family != "DTT" and is_TSD(dff, i):
            return i
    return np.nan

def get_TIR(x):
    s = x["seq"][200:-200]
    l_TIR = x["l_TIR"]
    return (s[0:l_TIR],s[-l_TIR:])

def get_TSD(x):
    s = x["seq"]
    l_TSD = x["l_TSD"]

    s1tsd = s[200 - l_TSD:200]
    last200 = s[-200:]
    s2tsd = last200[0:l_TSD]
    set1, set2 = sliding_window(s1tsd, s2tsd, l_TSD)
    tsd_dffset = get_difference(set1, set2)
    for i in tsd_dffset:
        if tsd_dffset[i] < l_TSD * 0.2:
            seq1 = set1[int(i.split(":")[0])]
            seq2 = set2[int(i.split(":")[1])]
            return (seq1, seq2)
    return np.nan


def TIR_TSD_percent(seq1, seq2):
    d = compare(seq1, seq2)
    l = len(seq1)
    p = (l - d) / l
    p = p * 100
    p = round(p, 2)
    return p


def process_result(df_in, module):
    df = df_in.copy()
    df["source"] = module
    df = df.loc[:, ["seqid", "source", "TIR_type", "sstart", "send", "TIR", "p_TIR", "TSD", "p_TSD", "len"]]
    df = df.rename(columns={"TIR_type": "type"})
    return df


def execute(df_in, module):
    df = df_in.copy()
    df["len"] = df.loc[:, "end"] - df.loc[:, "start"]
    df = df[df["len"] >= 450].reset_index(drop=True)
    df["l_TIR"] = df.swifter.progress_bar(True).apply(check_TIR, axis=1)
    df["l_TSD"] = df.swifter.progress_bar(True).apply(check_TSD, axis=1)

    df = df.dropna(ignore_index=True)
    df = df.astype({"l_TIR": int, "l_TSD": int})

    if(df.shape[0] == 0):
        return None

    df["TIR"] = df.swifter.progress_bar(True).apply(get_TIR, axis=1)
    df["p_TIR"] = df.swifter.progress_bar(True).apply(lambda x:
                                                      TIR_TSD_percent(x["TIR"][0],
                                                                      Seq(x["TIR"][1]).reverse_complement()), axis=1)
    df["TSD"] = df.swifter.progress_bar(True).apply(get_TSD, axis=1)
    df["p_TSD"] = df.swifter.progress_bar(True).apply(lambda x: TIR_TSD_percent(x["TSD"][0], x["TSD"][1]), axis=1)
    # both["rslt"] = both.swifter.progress_bar(True).apply(lambda x: x["seq"][200:-200], axis=1)
    df["len"] = df["len"] - 400
    return process_result(df, module)
