import math
import os
import shutil
import subprocess
import pandas as pd
import swifter  # ATTENTION: DO NOT REMOVE "swifter" EVEN IF IDE SHOWS IT IS NOT USED!
import multiprocessing as mp
from Bio import SeqIO

import prog_const


def prepare_fasta(genome_file, genome_name, GRF_mode, drop_seq_len):
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
    filtered_genome_file_name = f"{genome_name}{prog_const.spliter}filtered.fa"

    if GRF_mode == "boost":
        records_split_file_name = []
        os.makedirs("GRFmite_mp", exist_ok=True)
        os.chdir("./GRFmite_mp")
        for record in SeqIO.parse(genome_file, "fasta"):
            if len(record.seq) > drop_seq_len:
                file_name = f"{record.id}.fa"
                SeqIO.write(record, file_name, "fasta")
                records_split_file_name.append(file_name)
        os.chdir("../")
        return records_split_file_name, filtered_genome_file_name

    if GRF_mode == "mix":
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
        SeqIO.write(records_long, filtered_genome_file_name, "fasta")
        records_split_file_name = [f"{record.id}.fa" for record in records_short]
        return records_split_file_name, filtered_genome_file_name



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

    # mode = "native"
    records = []
    for record in SeqIO.parse(genome_file, "fasta"):
        if len(record.seq) > drop_seq_len:
            records.append(record)
    SeqIO.write(records, filtered_genome_file_name, "fasta")
    return None, filtered_genome_file_name


def GRF(GRF_path, file, cpu_cores, TIR_length):
    grf_bin_path = os.path.join(GRF_path, "grf-main")
    grf = (f"\"{grf_bin_path}\" -i \"{file}\" -o \"{file}_GRFmite\" -c 1 -t {int(cpu_cores)} -p 20 "
           f"--min_space 10 --max_space {int(TIR_length)} --max_indel 0 --min_tr 10 "
           f"--min_spacer_len 10 --max_spacer_len {int(TIR_length)}")
    shell_filter = r" | grep -vE 'start:|end:|print:|calculate|^$'"

    # TODO debug only
    # print(grf)
    # print(shell_filter)

    # add_three_spaces = r" | sed 's/^/   /'"
    # subprocess.Popen(grf + shell_filter + add_three_spaces, shell=True).wait()
    subprocess.Popen(grf + shell_filter, shell=True).wait()


def GRF_mp(GRF_path, fasta_files_list, TIR_length, num_GRF_instance, num_threads_per_GRF_instance):
    # To ensure one of the processes keep get <num_extra_threads> extra threads all the time requires
    # interprocess communication (IPC), which needs unnecessary complicated implementation for current task
    # To make the task simple, we just directly ensure every process get <num_extra_threads> extra threads
    # The worst case is having num_process * (num_process - 1) extra threads in total
    # which is generally acceptable in our task
    # mp_args_list = [(grfp, file, length, num_threads_per_GRF_instance + num_extra_threads) for file in fasta_files_list]
    mp_args_list = [(GRF_path, os.path.join("GRFmite_mp", file), num_threads_per_GRF_instance, TIR_length)
                    for file in fasta_files_list]
    # os.chdir("./GRFmite_mp")
    with mp.Pool(num_GRF_instance) as pool:
        pool.starmap(GRF, mp_args_list)
    # os.chdir("../")


def run_GRF_native(filtered_genome_file_name, GRF_path, cpu_cores, TIR_length):
    GRF(GRF_path, filtered_genome_file_name, cpu_cores, TIR_length)
    subprocess.Popen(["unlink", filtered_genome_file_name])


def run_GRF_mix(fasta_files_list, filtered_genome_file_name, GRF_path, cpu_cores, TIR_length):
    process_long_seq = mp.Process(target=GRF, args=(GRF_path, filtered_genome_file_name, cpu_cores, TIR_length))
    process_short_seq = mp.Process(target=GRF_mp, args=(GRF_path, fasta_files_list,
                                                        TIR_length, prog_const.mix_short_seq_process_num, 1))

    process_long_seq.start()
    process_short_seq.start()

    process_long_seq.join()
    process_short_seq.join()

    subprocess.Popen(["unlink", filtered_genome_file_name])
    subprocess.Popen(["rm", "-rf", "GRFmite_mp"])


def run_GRF_boost(fasta_files_list, GRF_result_dir_name, GRF_path, cpu_cores, TIR_length):
    num_processes, num_threads_per_process, num_extra_threads = cpu_cores_allocation_GRF_boost(cpu_cores, "cpu_bound")
    # print(num_processes, num_threads_per_process, num_extra_threads) # TODO debug only
    GRF_mp(GRF_path, fasta_files_list, TIR_length, num_processes, num_threads_per_process)
    collect_results(fasta_files_list, GRF_result_dir_name)
    subprocess.Popen(["rm", "-rf", "GRFmite_mp"])


def cpu_cores_allocation_GRF_boost(cpu_cores, job_bound_type="cpu_bound"):
    if job_bound_type == "cpu_bound":
        num_threads_total = cpu_cores + 1  # CPU Bound: N+1
    elif job_bound_type == "io_bound":
        num_threads_total = 2 * cpu_cores  # I/O Bound: 2N
    else:
        num_threads_total = cpu_cores

    num_processes = int(math.sqrt(num_threads_total))
    num_threads_per_process = int(num_threads_total / num_processes) * 4
    # num_extra_threads = num_threads_per_process - num_processes * num_threads_per_process
    # num_extra_threads = 0 if num_extra_threads < 0 else num_extra_threads
    num_extra_threads = 0
    return num_processes, num_threads_per_process, num_extra_threads


def collect_results(fasta_files_list, GRF_result_dir_name):
    GRF_result_files_list = [f"{file}_GRFmite" for file in fasta_files_list]
    os.makedirs(GRF_result_dir_name, exist_ok=True)
    with open(os.path.join(GRF_result_dir_name, "candidate.fasta"), "wb") as des:
        for f in GRF_result_files_list:
            try:
                with open(os.path.join("GRFmite_mp", f, "candidate.fasta"), "rb") as src:
                    shutil.copyfileobj(src, des)
            except FileNotFoundError:
                continue


def get_GRF_result_df(GRF_result_dir_name):
    GRF_result_file_path = os.path.join(GRF_result_dir_name, "candidate.fasta")

    # df_data_dict = {"id": [rec.id for rec in SeqIO.parse(GRF_result_file_path, "fasta")],
    # df_data_dict = {"id": list(SeqIO.index(GRF_result_file_path, "fasta")),
    #                 "seq": [str(rec.seq) for rec in SeqIO.parse(GRF_result_file_path, "fasta")],
    #                 "len": [len(rec) for rec in SeqIO.parse(GRF_result_file_path, "fasta")]}
    df_data_dict = [{"id": rec.id, "seq": str(rec.seq), "len": len(rec)}
                    for rec in SeqIO.parse(GRF_result_file_path, "fasta")]

    df_type = {"len": int}

    subprocess.Popen(["rm", "-rf", GRF_result_dir_name])
    return pd.DataFrame(df_data_dict, columns=["id", "seq", "len"]).astype(df_type)


def execute(TIRLearner_instance):
    genome_file = TIRLearner_instance.genome_file_path
    genome_name = TIRLearner_instance.genome_name
    cpu_cores = TIRLearner_instance.cpu_cores
    GRF_path = TIRLearner_instance.GRF_path
    TIR_length = TIRLearner_instance.TIR_length
    GRF_mode = TIRLearner_instance.GRF_mode

    drop_seq_len = int(TIR_length) + 500
    records_split_file_name, filtered_genome_file_name = prepare_fasta(genome_file, genome_name, GRF_mode, drop_seq_len)
    GRF_result_dir_name = f"{filtered_genome_file_name}_GRFmite"

    # print(os.listdir(os.getcwd()))
    # for i in os.listdir(os.getcwd()):
    #     shutil.copyfile(os.path.join(os.getcwd(), i), os.path.join(TIRLearner_instance.output_dir, i))

    print("  Step 1/2: Executing GRF\n")
    if GRF_mode == "mix":
        run_GRF_mix(records_split_file_name, filtered_genome_file_name, GRF_path, cpu_cores, TIR_length)
    elif GRF_mode == "boost":
        run_GRF_boost(records_split_file_name, GRF_result_dir_name, GRF_path, cpu_cores, TIR_length)
    else:
        run_GRF_native(filtered_genome_file_name, GRF_path, cpu_cores, TIR_length)
    print()

    print("  Step 2/2: Getting GRF result")
    return get_GRF_result_df(GRF_result_dir_name)
