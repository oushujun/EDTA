import math
import os
import shutil
import subprocess
import multiprocessing as mp
from Bio import SeqIO

import prog_const
spliter = prog_const.spliter


def prepare_fasta(genome_file, genome_name, mode, genome_num, drop_seq_len, short_seq_len=2000,
                  general_split_num_threshold=1, strict_split_percent_threshold=0.2):
    """This method performs two tasks:

      - Filters the input FASTA file by sequence length, removing sequences shorter than a specified filter length
        (filter_len), the filtered FASTA data will be used in the next step where we prepare specific fasta file(s)
        according to the mode;

      - Prepare specific fasta file(s) according to the mode:

        * native mode: the simplest way to execute GRF. One FASTA file is generated from the filtered FASTA data and
          served as the input of one GRF instance that will use the built-in multiprocessing function utilizing all the
          `t` available CPUs. This enables the first-layer multiprocessing.

        * boost mode: the most efficient way to execute GRF. Splits the filtered FASTA data into multiple FASTA files
          if the number of sequences exceeds a defined threshold (mp_threshold), otherwise redirect to the native mode.
          Each FASTA file contains only one sequence, every FASTA file will have its own GRF instance and be taken as
          the input, every GRF instance will use the built-in multiprocessing function utilizing all the `t` available
          CPUs. This enables the second-layer multiprocessing.

          By utilizing both first-layer multiprocessing (using the built-in multiprocessing function of GRF) and
          second-layer multiprocessing (creating a GRF instance for each processed FASTA file), we can optimize
          CPU usage for maximum efficiency. Gigantic amount (sequence number `genome_num` * available CPU number `t`)
          of threads will be created.

        * strict mode: the most strict way which tries to maximize efficiency under specific limit while executing GRF.
          The number of threads will be no bigger than the available CPU number `t`, ensuring one thread corresponds to
          one CPU (as much as possible, exclusions may apply). This mode is introduced to comply the strict CPU resource
          limitations on some HPC system to avoid penalties.

          Two pipelines of GRF execution will be created, each of them will have access to half of the available CPUs
          (1/2 * `t`). One pipeline will create one GRF instance and use the built-in multiprocessing function utilizing
          half of the available CPUs (determined by the available CPU number `t`) ...
          TODO documentation needs completion
    """
    if (mode == "native" or
            genome_num <= general_split_num_threshold or
            (mode == "strict" and drop_seq_len >= short_seq_len)):
        if mode != "native":
            print("   Number of sequences insufficient, " + mode + " mode unavailable, redirect to native mode")
        records = []
        for record in SeqIO.parse(genome_file, "fasta"):
            if len(record.seq) > drop_seq_len:
                records.append(record)
        SeqIO.write(records, f"{genome_name}_filtered.fa", "fasta")
        return None

    if mode == "strict":
        records_long = []
        records_short = []
        os.makedirs("GRFmite_mp", exist_ok=True)
        os.chdir("./GRFmite_mp")
        for record in SeqIO.parse(genome_file, "fasta"):
            record_len = len(record.seq)
            if record_len > drop_seq_len:
                if record_len > short_seq_len:
                    records_long.append(record)
                else:
                    records_short.append(record)
                    SeqIO.write(record, f"{record.id}.fa", "fasta")
        os.chdir("../")
        if len(records_short)/genome_num < strict_split_percent_threshold:
            subprocess.Popen(["rm", "-rf", "GRFmite_mp"])
            SeqIO.write(records_long + records_short, f"{genome_name}_filtered.fa", "fasta")
            print("   Percentage of short sequences insufficient, strict mode unavailable, redirect to native mode")
            return None

        SeqIO.write(records_long, f"{genome_name}_filtered.fa", "fasta")
        records_split_file_name = [f"{record.id}.fa" for record in records_short]
        return records_split_file_name

    records_split_file_name = []
    os.makedirs("GRFmite_mp", exist_ok=True)
    os.chdir("./GRFmite_mp")
    for record in SeqIO.parse(genome_file, "fasta"):
        if len(record.seq) > drop_seq_len:
            file_name = f"{record.id}.fa"
            SeqIO.write(record, file_name, "fasta")
            records_split_file_name.append(file_name)
    os.chdir("../")
    return records_split_file_name


def GRF(grfp, file, t, length):
    grf = (f"\"{grfp}/grf-main\" -i \"{file}\" -o \"{file}_GRFmite\" -c 1 -t {int(t)} -p 20 --min_space 10 "
           f"--max_space {int(length)} --max_indel 0 --min_tr 10 --min_spacer_len 10 --max_spacer_len {int(length)}")
    filter = r" | grep -vE 'start:|end:|print:|calculate|^$'"
    # add_three_spaces = r" | sed 's/^/   /'"
    # subprocess.Popen(grf + filter + add_three_spaces, shell=True).wait()
    subprocess.Popen(grf + filter, shell=True).wait()

def GRF_mp(grfp, fasta_files_list, t, length):
    mp_args_list = [(grfp, file, t, length) for file in fasta_files_list]
    os.chdir("./GRFmite_mp")
    with mp.Pool(int(t)) as pool:
        pool.starmap(GRF, mp_args_list)
    os.chdir("../")


def collect_results(genome_name, fasta_files_list):
    GRFmite_files_list = [f"{file}_GRFmite" for file in fasta_files_list]
    GRFmite_dir_name = f"{genome_name}_filtered.fa_GRFmite"
    os.makedirs(GRFmite_dir_name, exist_ok=True)
    with open(os.path.join(GRFmite_dir_name, "candidate.fasta"), "wb") as des:
        for f in GRFmite_files_list:
            try:
                with open(os.path.join("GRFmite_mp", f, "candidate.fasta"), "rb") as src:
                    shutil.copyfileobj(src, des)
            except FileNotFoundError:
                continue


def run_GRF_native(genome_name, grfp, t, length):
    GRF(grfp, f"{genome_name}_filtered.fa", t, length)
    subprocess.Popen(["unlink", f"{genome_name}_filtered.fa"])


def run_GRF_boost(fasta_files_list, genome_name, grfp, t, length):
        GRF_mp(grfp, fasta_files_list, t, length)
        collect_results(genome_name, fasta_files_list)
        subprocess.Popen(["rm", "-rf", "GRFmite_mp"])


# def GRF_strict_CPU_allocation(t):
#     short_seq_pipe_t = int(t/2)
#     long_seq_pipe_t = t - short_seq_pipe_t
#     short_seq_each_GRF_t = int(math.sqrt(short_seq_pipe_t))


def GRF_strict_CPU_allocation(t):
    short_seq_pipe_t = int(t/2)
    long_seq_pipe_t = t - short_seq_pipe_t




def run_GRF_strict(fasta_files_list, genome_name, grfp, t, length):
    process_long_seq = mp.Process(target=GRF, args=(grfp, f"{genome_name}_filtered.fa", t, length))
    process_short_seq = mp.Process(target=GRF_mp, args=(grfp, fasta_files_list, t, length))

    process_long_seq.start()
    process_short_seq.start()

    process_long_seq.join()
    process_short_seq.join()

    subprocess.Popen(["unlink", f"{genome_name}_filtered.fa"])
    subprocess.Popen(["rm", "-rf", "GRFmite_mp"])


def execute(TIRLearner_instance):
    genome_file = TIRLearner_instance.genome_file
    genome_name = TIRLearner_instance.genome_name
    t = TIRLearner_instance.cpu_cores
    grfp = TIRLearner_instance.GRF_path
    length = TIRLearner_instance.TIR_length
    mode = TIRLearner_instance.GRF_mode
    genome_num = TIRLearner_instance.genome_num

    drop_seq_len = int(length) + 500
    records_split_file_name = prepare_fasta(genome_file, genome_name, mode, genome_num, drop_seq_len)

    if records_split_file_name is None or mode == "native":
        run_GRF_native(genome_name, grfp, t, length)
    elif mode == "strict":
        # run_GRF_strict(records_split_file_name, genome_name, grfp, t, length)
        raise NotImplementedError
    else:
        run_GRF_boost(records_split_file_name, genome_name, grfp, t, length)
