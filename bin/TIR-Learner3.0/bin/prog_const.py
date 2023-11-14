import os

spliter = "-+-"
TIR_types = ("DTA", "DTC", "DTH", "DTM", "DTT")
program_root_dir = os.path.abspath(str(os.path.dirname(os.path.dirname(__file__))))
CNN_model_file = "cnn0912_tf_savedmodel"

print(program_root_dir)
