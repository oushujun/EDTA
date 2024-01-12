import os
import warnings

os.environ["TF_CPP_MIN_LOG_LEVEL"] = '3'  # mute all tensorflow info, warnings, and error msgs. #shujun
os.environ["KMP_WARNINGS"] = '0'  # mute all OpenMP warnings. #shujun
warnings.filterwarnings("ignore", category=FutureWarning)  # mute tensorflow warnings #shujun

# Use if True to suppress the PEP8: E402 warning
if True:  # noqa: E402
    import numpy as np
    import pandas as pd
    import swifter  # ATTENTION: DO NOT REMOVE "swifter" EVEN IF IDE SHOWS IT IS NOT USED!

    from sklearn.preprocessing import LabelEncoder
    # Attention: sklearn does not automatically import its subpackages
    import tensorflow as tf
    from keras.utils import to_categorical
    from keras.models import load_model

    import prog_const


def get_sequence_fragment(x, featureSize=200):
    seq = x["seq"]
    len_seq = len(seq)
    if len_seq >= featureSize * 2:
        return seq[0:featureSize] + seq[-featureSize:]
    s1 = seq[0:int(len_seq / 2)]
    s2 = seq[int(len_seq / 2):]
    n1 = "N" * (featureSize - len(s1))
    n2 = "N" * (featureSize - len(s2))
    s1 = s1 + n1
    s2 = n2 + s2
    return s1 + s2


def feature_encoding(df_in, flag_verbose):
    feature_int_encoder = LabelEncoder()
    voc = ["A", "C", "G", "T", "N"]
    num_classes = len(voc)
    feature_int_encoder.fit(voc)

    df = df_in.loc[:, ["id", "seq_frag"]].copy()
    print("  Step 2/8: Label Encoding - Transforming non-numerical labels to numerical labels")
    df["int_enc"] = df.swifter.progress_bar(flag_verbose).apply(
        lambda x: np.array(feature_int_encoder.transform(list(x["seq_frag"]))).reshape(-1, 1), axis=1)
    df = df.drop(columns="seq_frag")
    print("  Step 3/8: One-Hot Encoding - Converting class vectors to binary class matrices")
    df["feature"] = df.swifter.progress_bar(flag_verbose).apply(
        lambda x: to_categorical(x["int_enc"], num_classes=num_classes), axis=1)
    df = df.drop(columns="int_enc")

    # inputfeatures = np.array(input_features)
    # np.save(file + spliter + "features.npy", inputfeatures)
    return df


def predict(df_in, genome_file, path_to_model):
    model = load_model(path_to_model)
    pre_feature = df_in["feature"].to_numpy()
    df = df_in.drop(columns="feature")

    if pre_feature.shape[0] == 0:
        print("Info: " + genome_file + " has no candidate to be classified")
        return None

    l_class = ["DTA", "DTC", "DTH", "DTM", "DTT", "NonTIR"]
    target_int_encoder = LabelEncoder()
    target_int_encoder.fit(l_class)
    target_int_encoded = target_int_encoder.transform(l_class)
    d = dict(zip(target_int_encoded, l_class))

    print("  Step 5/8: Converting feature to tensor")
    with tf.device("/cpu:0"):
        pre_feature_tensor = tf.convert_to_tensor(np.stack(pre_feature), np.float32)

    predicted_labels = model.predict(pre_feature_tensor)
    df["percent"] = pd.Series(predicted_labels.max(axis=-1))
    y_classes = predicted_labels.argmax(axis=-1)
    df["TIR_type"] = pd.Series([d[i] for i in y_classes])
    return df


def postprocessing(df_in, flag_verbose):
    df = df_in.loc[:, ["id", "TIR_type"]]
    df = df[df["TIR_type"] != "NonTIR"].reset_index(drop=True)
    print("  Step 6/8: Retrieving sequence ID")
    df["seqid"] = df.swifter.progress_bar(flag_verbose).apply(lambda x: x["id"].split(":")[0], axis=1)
    print("  Step 7/8: Retrieving sequence starting coordinate")
    df["sstart"] = df.swifter.progress_bar(flag_verbose).apply(lambda x: int(x["id"].split(":")[1]), axis=1)
    print("  Step 8/8: Retrieving sequence ending coordinate")
    df["send"] = df.swifter.progress_bar(flag_verbose).apply(lambda x: int(x["id"].split(":")[2]), axis=1)
    df = df.loc[:, ["TIR_type", "id", "seqid", "sstart", "send"]]
    df = df.sort_values(["TIR_type", "seqid", "sstart", "send"], ignore_index=True)
    return df


def execute(TIRLearner_instance) -> pd.DataFrame:
    df = TIRLearner_instance["base"].copy()

    print("  Step 1/8: Getting sequence fragment for prediction")
    df["seq_frag"] = df.swifter.progress_bar(TIRLearner_instance.flag_verbose).apply(get_sequence_fragment, axis=1)
    df = df.drop(columns="seq")

    df = feature_encoding(df, TIRLearner_instance.flag_verbose)

    print("  Step 4/8: CNN prediction")
    df = predict(df, TIRLearner_instance.genome_file_path,
                 os.path.join(prog_const.program_root_dir_path, prog_const.CNN_model_dir_name))

    return postprocessing(df, TIRLearner_instance.flag_verbose)
