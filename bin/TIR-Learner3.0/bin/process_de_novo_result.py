import regex as re
import numpy as np
import pandas as pd
import swifter  # ATTENTION: DO NOT REMOVE "swifter" EVEN IF IDE SHOWS IT IS NOT USED!


def TA_repeats_check(df_in: pd.DataFrame, column: str = "seq", percent_threshold: float = 0.7) -> pd.DataFrame:
    df = df_in.copy()
    # df['TA_repeats_perc'] = (df["seq"].str.count('T') + df["seq"].str.count('A')) / df["seq"].str.len()
    df["TA_repeats_check"] = (df[column].str.count('T') + df[column].str.count('A') >=
                              df[column].str.len() * percent_threshold)
    return df[~df["TA_repeats_check"]].drop(columns="TA_repeats_check").reset_index(drop=True)


def check_N(s):
    n = s.count("N")
    if n > 0:
        return True
    return False


def check_N_per(s):
    n = s.count("N")
    p = n / len(s)
    if p >= 0.20:
        return True
    return False


def find_digits_sum(string):
    pattern = '(\d+)'
    l = re.findall(pattern, string)
    return sum([int(i) for i in l])


def TSD_check(x):
    TSD = x["TSD"]
    if ((len(TSD) > 6 or TSD == "TAA" or TSD == "TTA" or TSD == "TA" or x["seq"][0:4] == "CACT") or
            x["seq"][0:4] == "GTGA"):
        return True
    return False


def process_GRF_result(TIRLearner_instance):
    df_in = TIRLearner_instance.working_df_dict["GRF"]
    df = df_in[df_in["len"] >= 50].copy()

    if df.shape[0] == 0:
        print("NOTICE: No TIR found by GRF.")
        return None

    print("  Step 1/7: Getting TIR")
    df["TIR_len"] = df.swifter.progress_bar(TIRLearner_instance.flag_verbose).apply(
        lambda x: find_digits_sum(x["id"].split(":")[-2]), axis=1)
    # df["tirLen"] = df["tirLen"].astype(int)
    df["TIR"] = df.swifter.progress_bar(TIRLearner_instance.flag_verbose).apply(
        lambda x: x["seq"][0:x["TIR_len"]], axis=1)

    print("  Step 2/7: Checking TA repeats on sequence")
    # df["TA_repeats_seq_check"] = df.swifter.progress_bar(TIRLearner_instance.flag_verbose).apply(
    #     lambda x: np.nan if TA_repeats(x["seq"]) else False, axis=1)
    # df = df.dropna(ignore_index=True)
    df = TA_repeats_check(df)

    print("  Step 3/7: Checking percentage of N on sequence")
    df["check_N_per_seq_check"] = df.swifter.progress_bar(TIRLearner_instance.flag_verbose).apply(
        lambda x: np.nan if check_N_per(x["seq"]) else False, axis=1)
    df = df.dropna(ignore_index=True).copy()

    print("  Step 4/7: Checking TA repeats on TIR")
    # df["TA_repeats_tir_check"] = df.swifter.progress_bar(TIRLearner_instance.flag_verbose).apply(
    #     lambda x: np.nan if TA_repeats(x["tir"]) else False, axis=1)
    # df = df.dropna(ignore_index=True)
    df = TA_repeats_check(df, "TIR")

    print("  Step 5/7: Checking N existence on TIR")
    df["check_N_TIR_check"] = df.swifter.progress_bar(TIRLearner_instance.flag_verbose).apply(
        lambda x: np.nan if check_N(x["TIR"]) else False, axis=1)
    df = df.dropna(ignore_index=True).copy()

    print("  Step 6/7: Getting TSD")
    df["TSD"] = df.swifter.progress_bar(TIRLearner_instance.flag_verbose).apply(
        lambda x: x["id"].split(":")[-1], axis=1)

    print("  Step 7/7: Checking TSD")
    df["TSD_check"] = df.swifter.progress_bar(TIRLearner_instance.flag_verbose).apply(
        lambda x: True if TSD_check(x) else np.nan, axis=1)
    df = df.dropna(ignore_index=True).loc[:, ["id", "seq"]].copy()

    # print("  Step 9/9: Saving processed GRFmite")
    # df["id"] = df.swifter.progress_bar(TIRLearner_instance.flag_verbose).apply(lambda x: ">" + x["id"], axis=1)
    df["id"] = ">" + df["id"]
    # df.to_csv(os.path.join("../", f"{genome_name}{spliter}processedGRFmite.fa"), sep='\n', header=False, index=False)
    # df.to_csv(TIRLearner_instance.processed_de_novo_result_file, sep='\n', header=False, index=False)
    return df


def process_TIRvish_result(TIRLearner_instance):
    df_in = TIRLearner_instance["TIRvish"]
    df = df_in[df_in["end"] - df_in["start"] + 1 >= 50].copy()

    if df.shape[0] == 0:
        print("NOTICE: No TIR found by TIRvish.")
        return None

    print("  Step 1/5: Getting TIR")
    df["TIR1_start"] = df["TIR1_start"] - df["start"]
    df.loc[df["TIR1_start"] < 0, "TIR1_start"] = 0
    df["TIR1_end"] = df["TIR1_end"] - df["start"]

    df["TIR2_start"] = df["TIR2_start"] - df["start"]
    df.loc[df["TIR2_start"] < 0, "TIR2_start"] = 0
    df["TIR2_end"] = df["TIR2_end"] - df["start"]

    df["TIR1"] = df.swifter.progress_bar(TIRLearner_instance.flag_verbose).apply(
        lambda x: x["seq"][x["TIR1_start"]: x["TIR1_end"]], axis=1)
    df["TIR2"] = df.swifter.progress_bar(TIRLearner_instance.flag_verbose).apply(
        lambda x: x["seq"][x["TIR2_start"]: x["TIR2_end"]], axis=1)

    print("  Step 2/5: Checking TA repeats on sequence")
    df = TA_repeats_check(df)

    print("  Step 3/5: Checking percentage of N on sequence")
    df["check_N_per_seq_check"] = df.swifter.progress_bar(TIRLearner_instance.flag_verbose).apply(
        lambda x: np.nan if check_N_per(x["seq"]) else False, axis=1)
    df = df.dropna(ignore_index=True).copy()

    print("  Step 4/5: Checking TA repeats on TIR")
    df = TA_repeats_check(df, "TIR1")
    df = TA_repeats_check(df, "TIR2")

    print("  Step 5/5: Checking N existance on TIR")
    df["check_N_TIR1_check"] = df.swifter.progress_bar(TIRLearner_instance.flag_verbose).apply(
        lambda x: np.nan if check_N(x["TIR1"]) else False, axis=1)
    df = df.dropna(ignore_index=True).copy()
    df["check_N_TIR2_check"] = df.swifter.progress_bar(TIRLearner_instance.flag_verbose).apply(
        lambda x: np.nan if check_N(x["TIR2"]) else False, axis=1)
    df = df.dropna(ignore_index=True).loc[:, ["id", "seq"]].copy()

    # print("  Step 6/6: Saving processed TIRvish")
    # df.to_csv(TIRLearner_instance.processed_de_novo_result_file, sep='\n', header=False, index=False)
    return df


def combine_de_novo_result(TIRLearner_instance):
    df = pd.concat((TIRLearner_instance.get("TIRvish", None), TIRLearner_instance.get("GRF", None)), ignore_index=True)
    df.to_csv(TIRLearner_instance.processed_de_novo_result_file_name, sep='\n', header=False, index=False)
    del TIRLearner_instance["TIRvish"]
    del TIRLearner_instance["GRF"]
