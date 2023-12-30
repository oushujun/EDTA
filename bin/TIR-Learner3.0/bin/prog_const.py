import os

spliter = "-+-"
TIR_types = ("DTA", "DTC", "DTH", "DTM", "DTT")
program_root_dir = os.path.abspath(str(os.path.dirname(os.path.dirname(__file__))))
CNN_model_file = "cnn0912_tf_savedmodel"

short_seq_len = 2000
general_split_num_threshold = 5
mix_split_percent_threshold = 0.05
