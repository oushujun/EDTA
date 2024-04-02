# import os
# import subprocess
# import pandas as pd
# import swifter  # ATTENTION: DO NOT REMOVE "swifter" EVEN IF IDE SHOWS IT IS NOT USED!
# from get_fasta_sequence import get_fasta_pieces_SeqIO
#
# import prog_const

from const import *

from get_fasta_sequence import get_fasta_pieces_SeqIO


def split_sequence_evenly(seq_record, split_seq_len, overlap_seq_len):
    """
    Evenly splits a sequence into several segmented sequences, each segment has a maximum length of
    the length threshold split_seq_len.

    Parameters:
    - seq_record: A SeqRecord object representing the sequence.
    - split_seq_len: The length threshold for splitting, i.e. maximum allowed sequence length before splitting.
    - overlap_seq_len: An overlapping segment that covers the boundary of the preceding segment and
    the succeeding segment with the length of overlap_seq_len will be created. overlap_seq_len // 2 represents
    the length of the two parts that come from the preceding segment and the succeeding segment.
    Specifying a value <=1 will result in no overlapping segment being produced.

    Returns:
    A list of SeqRecord objects representing the original or split sequences.
    """
    records = []
    seq_len = len(seq_record)
    overlap_seg_half_len = overlap_seq_len // 2
    total_parts = seq_len // split_seq_len + (1 if seq_len % split_seq_len else 0)
    for i in range(total_parts):
        segment_start = i * split_seq_len
        segment_end = min(segment_start + split_seq_len, seq_len)
        segment_seq = seq_record.seq[segment_start:segment_end]
        part_id = f"{seq_record.id}_split_{i + 1}of{total_parts}"
        records.append(SeqRecord(segment_seq, id=part_id, description=""))
        if overlap_seg_half_len > 0 and segment_end != seq_len:
            overlap_seg = seq_record.seq[segment_end - overlap_seg_half_len:
                                         min(segment_end + overlap_seg_half_len, seq_len)]
            overlap_id = f"{seq_record.id}_split_{i + 1}.5of{total_parts}"
            records.append(SeqRecord(overlap_seg, id=overlap_id, description=""))
    return records


# def process_fasta(genome_file):
#     """
#     Write each sequence in the FASTA file into separate FASTA files and further split the sequence into segments if
#     when needed based on a length threshold split_seq_len.
#
#     Parameters:
#     - file_name: Path to the FASTA file.
#     - split_seq_len: Length threshold for splitting sequence.
#     """
#     split_fasta_files_list = []
#     for seq_record in SeqIO.parse(genome_file, "fasta"):
#         if len(seq_record) >= tirvish_split_seq_len:
#             segments = split_sequence_evenly(seq_record, tirvish_split_seq_len, tirvish_overlap_seq_len)
#             for segment in segments:
#                 file_name = f"{segment.id}.fasta"
#                 SeqIO.write(segment, file_name, "fasta")
#                 split_fasta_files_list.append(file_name)
#         else:
#             file_name = f"{seq_record.id}.fasta"
#             SeqIO.write(seq_record, file_name, "fasta")
#             split_fasta_files_list.append(file_name)
#     return split_fasta_files_list


def save_fasta_file(seq_record):
    fasta_file_name = f"{seq_record.id}.fasta"
    working_dir_name = f"{fasta_file_name}_{splited_fasta_tag}"
    os.makedirs(working_dir_name, exist_ok=True)
    fasta_file_path = os.path.join(working_dir_name, fasta_file_name)
    SeqIO.write(seq_record, fasta_file_path, "fasta")
    return fasta_file_path


def process_fasta(genome_file, split_seq_len, overlap_seq_len):
    # TODO: Universal process fasta function for both TIRvish and GRF?
    """
    Write each sequence in the FASTA file into separate FASTA files and further split the sequence into segments if
    when needed based on a length threshold split_seq_len.

    Parameters:
    - file_name: Path to the FASTA file.
    - split_seq_len: Length threshold for splitting sequence.
    """
    split_fasta_files_path_list = []
    for seq_record in SeqIO.parse(genome_file, "fasta"):
        if len(seq_record) >= TIRvish_split_seq_len:
            segments = split_sequence_evenly(seq_record, split_seq_len, overlap_seq_len)
            for segment in segments:
                split_fasta_files_path_list.append(save_fasta_file(segment))
        else:
            split_fasta_files_path_list.append(save_fasta_file(seq_record))
    return split_fasta_files_path_list


def retrieve_split_sequence_offset(segment_position, split_seq_len, overlap_seq_len):
    offset = 0
    if split_seq_len == 0:
        raise ValueError("When split happens, split_seq_len cannot be zero.")
    try:
        # Normal segment, segment_position is like: <x>of<y>, where x and y are integers
        segment_index = int(segment_position.split("of")[0])
        offset = (segment_index - 1) * split_seq_len
    except ValueError:
        # Overlap segment, segment_position is like: <x>.5of<y>, where x and y are integers
        segment_index = int(segment_position.split("of")[0][:-2])  # remove the .5
        offset = segment_index * split_seq_len - overlap_seq_len // 2
    return offset


def TIRvish(genome_file, genome_name, TIR_length, gt_path):
    # print(os.listdir())  # TODO only for debug
    gt_bin_path = os.path.join(gt_path, "gt")
    gt_index_file_name = genome_name + spliter + "gt_index"
    subprocess.Popen(
        [gt_bin_path, "suffixerator", "-db", genome_file, "-indexname", gt_index_file_name,
         "-tis", "-suf", "-lcp", "-des", "-ssp", "-sds", "-dna", "-mirrored"]).wait()

    TIRvish_result_gff3_file_name = f"{genome_name}{spliter}TIRvish.gff3"
    gt_tirvish = (f"\"{gt_bin_path}\" tirvish -index {gt_index_file_name} -seed 20 -mintirlen 10 -maxtirlen 1000 "
                  f"-mintirdist 10 -maxtirdist {str(TIR_length)} -similar 80 -mintsd 2 -maxtsd 11 "
                  f"-vic 13 -seqids \"yes\" > {TIRvish_result_gff3_file_name}")
    subprocess.Popen(gt_tirvish, shell=True).wait()
    subprocess.Popen(["find", ".", "-name", f"{gt_index_file_name}*", "-delete"])
    return TIRvish_result_gff3_file_name


def TIRvish_mp(genome_file_path, genome_name, TIR_length, gt_path):
    TIRvish_working_dir = os.path.dirname(genome_file_path)
    genome_file_name = os.path.basename(genome_file_path)
    os.chdir(TIRvish_working_dir)
    TIRvish_result_gff3_file_name = TIRvish(genome_file_name, genome_name, TIR_length, gt_path)
    os.chdir("../")
    return os.path.join(TIRvish_working_dir, TIRvish_result_gff3_file_name)


def get_TIRvish_result_df(TIRvish_result_gff3_file_path, flag_debug, split_seq_len=0, overlap_seq_len=0):
    df_data_dict = {"seqid": [], "start": [], "end": [], "TIR1_start": [], "TIR1_end": [], "TIR2_start": [],
                    "TIR2_end": [], "id": []}
    df_type = {"start": int, "end": int, "TIR1_start": int, "TIR1_end": int, "TIR2_start": int, "TIR2_end": int}

    if os.path.exists(TIRvish_result_gff3_file_path) and os.path.getsize(TIRvish_result_gff3_file_path) != 0:
        with open(TIRvish_result_gff3_file_path, 'r') as f:
            while (line := f.readline()) != "":
                while line.startswith('#'):
                    line = f.readline()
                TSD1 = list(map(int, f.readline().split('\t')[3:5]))
                record = f.readline().split('\t')
                TIR1 = list(map(int, f.readline().split('\t')[3:5]))
                TIR2 = list(map(int, f.readline().split('\t')[3:5]))
                TSD2 = list(map(int, f.readline().split('\t')[3:5]))
                f.readline()  # Jump the line "###"

                seqid = record[0]
                start = int(record[3])
                end = int(record[4])

                offset = 0
                # Offset for segment of sequence from splitting
                if "_split_" in seqid:
                    # seqid = <original_seqid>_split_<segment_position>
                    segment_position = seqid.split("_split_")[1]
                    seqid = seqid.split("_split_")[0]  # reassign original_seqid to seqid
                    offset = retrieve_split_sequence_offset(segment_position, split_seq_len, overlap_seq_len)

                # -1 at start coordinates are internal adjustments for get_fasta_sequence of TIR-Learner
                # Actual start coordinates are without the -1 adjustments, i.e. start + offset
                df_data_dict["seqid"].append(seqid)
                df_data_dict["start"].append(start - 1 + offset)
                df_data_dict["end"].append(end + offset)
                df_data_dict["TIR1_start"].append(TIR1[0] - 1 + offset)
                df_data_dict["TIR1_end"].append(TIR1[1] + offset)
                df_data_dict["TIR2_start"].append(TIR2[0] - 1 + offset)
                df_data_dict["TIR2_end"].append(TIR2[1] + offset)
                df_data_dict["id"].append(
                    f">{seqid}:{start + offset}:{end + offset}:"
                    f"tsd1_{TSD1[0] + offset}_{TSD1[1] + offset}_tsd2_{TSD2[0] + offset}_{TSD2[1] + offset}_"
                    f"tir1_{TIR1[0] + offset}_{TIR1[1] + offset}_tir2_{TIR2[0] + offset}_{TIR2[1] + offset}")

    if not flag_debug:
        subprocess.Popen(["unlink", TIRvish_result_gff3_file_path])
    return pd.DataFrame(df_data_dict).astype(df_type)


def run_TIRvish_native(genome_file, genome_name, TIR_length, flag_debug, gt_path):
    print("  Step 1/2: Executing TIRvish in native mode")
    TIRvish_result_gff3_file_name = TIRvish(genome_file, genome_name, TIR_length, gt_path)
    print("  Step 2/2: Getting TIRvish result")
    return get_TIRvish_result_df(TIRvish_result_gff3_file_name, flag_debug)


def run_TIRvish_py_para(genome_file, genome_name, TIR_length, cpu_cores, flag_debug, gt_path, fasta_files_path_list):
    os.makedirs(f"{splited_fasta_tag}_mp", exist_ok=True)
    os.chdir(f"./{splited_fasta_tag}_mp")

    print("  Step 1/3: Processing FASTA files")
    fasta_files_path_list.extend(process_fasta(genome_file, TIRvish_split_seq_len, TIRvish_overlap_seq_len))

    print("  Step 2/3: Executing TIRvish with python multiprocess")
    mp_args_list = [(file_path, genome_name, TIR_length, gt_path) for file_path in fasta_files_path_list]
    with mp.Pool(cpu_cores) as pool:
        TIRvish_result_gff3_file_path_list = pool.starmap(TIRvish_mp, mp_args_list)

    print("  Step 3/3: Getting TIRvish result")
    mp_args_list = [(file_path, flag_debug, TIRvish_split_seq_len, TIRvish_overlap_seq_len) for file_path in
                    TIRvish_result_gff3_file_path_list]
    with mp.Pool(cpu_cores * thread_core_ratio) as pool:
        TIRvish_result_df_list = pool.starmap(get_TIRvish_result_df, mp_args_list)

    os.chdir("../")
    return pd.concat(TIRvish_result_df_list).reset_index(drop=True)


def execute(TIRLearner_instance):
    genome_file = TIRLearner_instance.genome_file_path
    genome_name = TIRLearner_instance.genome_name
    TIR_length = TIRLearner_instance.TIR_length
    cpu_cores = TIRLearner_instance.cpu_cores
    para_mode = TIRLearner_instance.para_mode
    # TODO add GNU Parallel support
    flag_verbose = TIRLearner_instance.flag_verbose
    flag_debug = TIRLearner_instance.flag_debug
    gt_path = TIRLearner_instance.gt_path
    additional_args = TIRLearner_instance.additional_args
    fasta_files_path_list = TIRLearner_instance.split_fasta_files_path_list

    if NO_PARALLEL in additional_args:
        df = run_TIRvish_native(genome_file, genome_name, TIR_length, flag_debug, gt_path)
    elif para_mode == "gnup":
        raise NotImplementedError()
    else:
        df = run_TIRvish_py_para(genome_file, genome_name, TIR_length, cpu_cores, flag_debug, gt_path,
                                 fasta_files_path_list)

    return get_fasta_pieces_SeqIO(genome_file, df, cpu_cores, flag_verbose)
