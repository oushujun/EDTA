import os

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


def process_additional_args(additional_args: list) -> tuple:
    processed_additional_args = tuple(i for i in
                                      tuple(map(additional_args_mapping_dict.get, additional_args)) if i is not None)
    if SKIP_TIRVISH in processed_additional_args and SKIP_GRF in processed_additional_args:
        raise SystemExit("ERROR: \"skip_tirvish\" and \"skip_grf\" cannot be specified at the same time!")
    return processed_additional_args
