import os
import warnings

os.environ["TF_CPP_MIN_LOG_LEVEL"] = '3'  # mute all tensorflow info, warnings, and error msgs. #shujun
os.environ["KMP_WARNINGS"] = '0'  # mute all OpenMP warnings. #shujun
warnings.filterwarnings("ignore", category=FutureWarning)  # mute tensorflow warnings and pyarrow warning
warnings.filterwarnings("ignore", category=UserWarning)  # mute h5py warning

# Use if True to suppress the PEP8: E402 warning
if True:  # noqa: E402
    import datetime
    import json
    import math
    import multiprocessing as mp
    import regex as re
    import shutil
    import subprocess
    import tempfile
    import time

    import numpy as np
    import pandas as pd
    import swifter

    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    from sklearn.preprocessing import LabelEncoder
    # Attention: sklearn does not automatically import its subpackages
    import tensorflow as tf
    from tensorflow.python.framework.errors_impl import InternalError
    from keras.utils import to_categorical
    from keras.models import load_model


# Acceptable additional args
CHECKPOINT_OFF = "CHECKPOINT_OFF"
NO_PARALLEL = "NO_PARALLEL"
SKIP_TIRVISH = "SKIP_TIRVISH"
SKIP_GRF = "SKIP_GRF"

spliter = "-+-"
TIR_types = ("DTA", "DTC", "DTH", "DTM", "DTT")

CNN_model_dir_name = "cnn0912_tf_savedmodel"
sandbox_dir_name = "[DONT_ALTER]TIR-Learner_sandbox"
splited_fasta_tag = "SplitedFasta"

program_root_dir_path = os.path.abspath(str(os.path.dirname(os.path.dirname(__file__))))

ref_lib_dir_name = "RefLib"
ref_lib_available_species = ("rice", "maize")
ref_lib_file_dict = {species: [f"{species}_{TIR_type}_RefLib" for TIR_type in TIR_types]
                     for species in ref_lib_available_species}
ref_lib_dir_path = os.path.join(program_root_dir_path, ref_lib_dir_name)

TIRvish_split_seq_len = 5 * (10 ** 6)  # 5 mb
TIRvish_overlap_seq_len = 50 * (10 ** 3)  # 50 kb

# TODO only for debug
# TIRvish_split_seq_len = 10000
# TIRvish_overlap_seq_len = 10000


short_seq_len = 2000
general_split_num_threshold = 5
mix_split_percent_threshold = 0.05
mix_short_seq_process_num = 2

thread_core_ratio = 2


def process_additional_args(additional_args: list) -> tuple:
    if additional_args == [""]:
        return tuple()
    processed_additional_args = tuple(map(str.upper, additional_args))
    if (SKIP_TIRVISH in processed_additional_args) and (SKIP_GRF in processed_additional_args):
        raise SystemExit("ERROR: \"skip_tirvish\" and \"skip_grf\" cannot be specified at the same time!")
    return processed_additional_args
