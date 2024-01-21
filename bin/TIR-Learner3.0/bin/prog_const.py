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

    import numpy as np
    import pandas as pd
    import swifter

    from Bio import SeqIO
    from Bio.Seq import Seq

    from sklearn.preprocessing import LabelEncoder
    # Attention: sklearn does not automatically import its subpackages
    import tensorflow as tf
    from tensorflow.python.framework.errors_impl import InternalError
    from keras.utils import to_categorical
    from keras.models import load_model


# Acceptable additional args
FORCE_GRF_MODE = "FORCE_GRF_MODE"
CHECKPOINT_OFF = "CHECKPOINT_OFF"
SKIP_TIRVISH = "SKIP_TIRVISH"
SKIP_GRF = "SKIP_GRF"

spliter = "-+-"
TIR_types = ("DTA", "DTC", "DTH", "DTM", "DTT")
additional_args_mapping_dict = {"force_grf_mode": FORCE_GRF_MODE,
                                "checkpoint_off": CHECKPOINT_OFF,
                                "skip_tirvish": SKIP_TIRVISH,
                                "skip_grf": SKIP_GRF}

CNN_model_dir_name = "cnn0912_tf_savedmodel"
sandbox_dir_name = "[DO_NOT_ALTER]_TIR-Learner_sandbox_directory"

program_root_dir_path = os.path.abspath(str(os.path.dirname(os.path.dirname(__file__))))

ref_lib_dir_name = "RefLib"
ref_lib_available_species = ("rice", "maize")
ref_lib_file_dict = {species: [f"{species}_{TIR_type}_RefLib" for TIR_type in TIR_types]
                     for species in ref_lib_available_species}
ref_lib_dir_path = os.path.join(program_root_dir_path, ref_lib_dir_name)

short_seq_len = 2000
general_split_num_threshold = 5
mix_split_percent_threshold = 0.05
mix_short_seq_process_num = 2

thread_core_ratio = 2


def process_additional_args(additional_args: list) -> tuple:
    processed_additional_args = tuple(i for i in
                                      tuple(map(additional_args_mapping_dict.get, additional_args)) if i is not None)
    if SKIP_TIRVISH in processed_additional_args and SKIP_GRF in processed_additional_args:
        raise SystemExit("ERROR: \"skip_tirvish\" and \"skip_grf\" cannot be specified at the same time!")
    return processed_additional_args
