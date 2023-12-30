import math
import os
import shutil
import subprocess
import multiprocessing as mp
from Bio import SeqIO

# from typing import TYPE_CHECKING
#
# if TYPE_CHECKING:
#     from main import TIRLearner

import prog_const


def prepare_fasta(genome_file, genome_name, mode, drop_seq_len):
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

        * mix mode: the most strict way which tries to maximize efficiency under specific limit while executing GRF.
          The number of threads will be no bigger than the available CPU number `t`, ensuring one thread corresponds to
          one CPU (as much as possible, exclusions may apply). This mode is introduced to comply the strict CPU resource
          limitations on some HPC system to avoid penalties.

          Two pipelines of GRF execution will be created, each of them will have access to half of the available CPUs
          (1/2 * `t`). One pipeline will create one GRF instance and use the built-in multiprocessing function utilizing
          half of the available CPUs (determined by the available CPU number `t`) ...
          TODO documentation needs complete revision!
    """
    if mode == "native":
        records = []
        for record in SeqIO.parse(genome_file, "fasta"):
            if len(record.seq) > drop_seq_len:
                records.append(record)
        SeqIO.write(records, f"{genome_name}_filtered.fa", "fasta")
        return None

    if mode == "mix":
        records_long = []
        records_short = []
        os.makedirs("GRFmite_mp", exist_ok=True)
        os.chdir("./GRFmite_mp")
        for record in SeqIO.parse(genome_file, "fasta"):
            record_len = len(record.seq)
            if record_len > drop_seq_len:
                if record_len > prog_const.short_seq_len:
                    records_long.append(record)
                else:
                    records_short.append(record)
                    SeqIO.write(record, f"{record.id}.fa", "fasta")
        os.chdir("../")
        SeqIO.write(records_long, f"{genome_name}_filtered.fa", "fasta")
        records_split_file_name = [f"{record.id}.fa" for record in records_short]
        return records_split_file_name

        # if cpu_cores >= mix_short_seq_process_num / mix_split_percent_threshold:
        #     # Mix mode will only available when cpu_cores >= 10,
        #     # which means we can allocate 20% of the cpu at most to process the short sequences
        #     print(f"   Number of available CPU cores insufficient "
        #           f"(should be at least {mix_short_seq_process_num / mix_split_percent_threshold}), "
        #           f"{mode} mode unavailable, redirect to native mode.")
        #     print("   Boost mode may be applicable, if you wish to execute in boost mode, "
        #           "please rerun the program and specify \"-m boost\"")
        #     subprocess.Popen(["rm", "-rf", "GRFmite_mp"])
        #     SeqIO.write(records_long + records_short, f"{genome_name}_filtered.fa", "fasta")
        #     return None
        #
        # if len(records_short) / genome_num < mix_split_percent_threshold:
        #     print(f"   Percentage of short sequences insufficient "
        #           f"(should be at least {mix_split_percent_threshold * 100}%), "
        #           f"{mode} mode unneeded, redirect to native mode")
        #     # print("   Boost mode may be applicable, if you wish to execute in boost mode, "
        #     #       "please rerun the program and specify \"-m boost\"")
        #     subprocess.Popen(["rm", "-rf", "GRFmite_mp"])
        #     SeqIO.write(records_long + records_short, f"{genome_name}_filtered.fa", "fasta")
        #     return None
        #
        # if len(records_short)/genome_num > (1 - mix_split_percent_threshold):
        #     print(f"   WARNING: Two many short sequences exist (over {(1 - mix_split_percent_threshold) * 100}%), "
        #           f"you will very likely not getting useful result!")

    # mode = "boost"
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
    shell_filter = r" | grep -vE 'start:|end:|print:|calculate|^$'"

    # TODO debug only
    # print(grf)
    # print(shell_filter)

    # add_three_spaces = r" | sed 's/^/   /'"
    # subprocess.Popen(grf + shell_filter + add_three_spaces, shell=True).wait()
    subprocess.Popen(grf + shell_filter, shell=True).wait()


def GRF_mp(grfp, fasta_files_list, length, num_GRF_instance, num_threads_per_GRF_instance, num_extra_threads=0):
    # To ensure one of the processes keep get <num_extra_threads> extra threads all the time requires
    # interprocess communication (IPC), which needs unnecessary complicated implementation for current task
    # To make the task simple, we just directly ensure every process get <num_extra_threads> extra threads
    # The worst case is having num_process * (num_process - 1) extra threads in total
    # which is generally acceptable in our task
    # mp_args_list = [(grfp, file, length, num_threads_per_GRF_instance + num_extra_threads) for file in fasta_files_list]
    mp_args_list = [(grfp, os.path.join("GRFmite_mp", file), num_threads_per_GRF_instance + num_extra_threads, length)
                    for file in fasta_files_list]
    # os.chdir("./GRFmite_mp")
    with mp.Pool(num_GRF_instance) as pool:
        pool.starmap(GRF, mp_args_list)
    # os.chdir("../")


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


def cpu_cores_allocation_GRF_boost(cpu_cores, job_bound_type="cpu_bound"):
    if job_bound_type == "cpu_bound":
        num_threads_total = cpu_cores + 1  # CPU Bound: N+1
    elif job_bound_type == "io_bound":
        num_threads_total = 2 * cpu_cores  # I/O Bound: 2N
    else:
        num_threads_total = cpu_cores

    num_processes = int(math.sqrt(num_threads_total))
    num_threads_per_process = int(num_threads_total / num_processes)
    num_extra_threads = num_threads_per_process - num_processes * num_threads_per_process
    num_extra_threads = 0 if num_extra_threads < 0 else num_extra_threads
    return num_processes, num_threads_per_process, num_extra_threads


def run_GRF_boost(fasta_files_list, genome_name, grfp, cpu_cores, length):
    num_processes, num_threads_per_process, num_extra_threads = cpu_cores_allocation_GRF_boost(cpu_cores, "cpu_bound")
    # print(num_processes, num_threads_per_process, num_extra_threads) # TODO debug only
    GRF_mp(grfp, fasta_files_list, length, num_processes, num_threads_per_process, num_extra_threads)
    collect_results(genome_name, fasta_files_list)
    subprocess.Popen(["rm", "-rf", "GRFmite_mp"])


def run_GRF_mix(fasta_files_list, genome_name, grfp, cpu_cores, length, mix_short_seq_process_num=2):
    process_long_seq = mp.Process(target=GRF, args=(grfp, f"{genome_name}_filtered.fa", cpu_cores, length))
    process_short_seq = mp.Process(target=GRF_mp, args=(grfp, fasta_files_list, length, mix_short_seq_process_num, 1))

    process_long_seq.start()
    process_short_seq.start()

    process_long_seq.join()
    process_short_seq.join()

    subprocess.Popen(["unlink", f"{genome_name}_filtered.fa"])
    subprocess.Popen(["rm", "-rf", "GRFmite_mp"])


def execute(TIRLearner_instance):
    genome_file = TIRLearner_instance.genome_file
    genome_name = TIRLearner_instance.genome_name
    cpu_cores = TIRLearner_instance.cpu_cores
    grfp = TIRLearner_instance.GRF_path
    length = TIRLearner_instance.TIR_length
    mode = TIRLearner_instance.GRF_mode

    drop_seq_len = int(length) + 500
    records_split_file_name = prepare_fasta(genome_file, genome_name, mode, drop_seq_len)

    # print(os.listdir(os.getcwd()))
    # for i in os.listdir(os.getcwd()):
    #     shutil.copyfile(os.path.join(os.getcwd(), i), os.path.join(TIRLearner_instance.output_dir, i))

    if records_split_file_name is None or mode == "native":
        run_GRF_native(genome_name, grfp, cpu_cores, length)
    elif mode == "mix":
        run_GRF_mix(records_split_file_name, genome_name, grfp, cpu_cores, length)
    else:
        run_GRF_boost(records_split_file_name, genome_name, grfp, cpu_cores, length)
