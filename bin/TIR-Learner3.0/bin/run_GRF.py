import os
import shutil
import subprocess
import multiprocessing as mp
from Bio import SeqIO

spliter = "-+-"


def filter_fasta_by_length_split_when_mp_needed(genome_file, genome_name, filter_len, mp_threshold=1):
    num_seq = len([True for line in open(genome_file) if line.startswith(">")])
    records = []
    if num_seq > mp_threshold:
        os.makedirs("GRFmite_mp", exist_ok=True)
        os.chdir("./GRFmite_mp")
        for record in SeqIO.parse(genome_file, "fasta"):
            if len(record.seq) > filter_len:
                file_name = f"{record.id}.fa"
                SeqIO.write(record, file_name, "fasta")
                records.append(file_name)
        os.chdir("../")
    else:
        for record in SeqIO.parse(genome_file, "fasta"):
            if len(record.seq) > filter_len:
                records.append(record)
        SeqIO.write(records, f"{genome_name}_filtered.fa", "fasta")
        records = None
    return records


def GRF(grfp, file, t, length):
    grf = (f"\"{grfp}/grf-main\" -i \"{file}\" -o \"{file}_GRFmite\" -c 1 -t {int(t)} -p 20 --min_space 10 "
           f"--max_space {int(length)} --max_indel 0 --min_tr 10 --min_spacer_len 10 --max_spacer_len {int(length)}")
    filter = r" | grep -vE 'start:|end:|print:|calculate|^$'"
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


def execute(args):
    genome_file = args[0]
    genome_name = args[1]
    t = args[3]
    grfp = args[5]
    length = args[6]

    filter_len = int(length) + 500
    fasta_files_list = filter_fasta_by_length_split_when_mp_needed(genome_file, genome_name, filter_len)

    print("Running GRF")
    if fasta_files_list is None:
        GRF(grfp, f"{genome_name}_filtered.fa", t, length)
        subprocess.Popen(["unlink", f"{genome_name}_filtered.fa"])
    else:
        GRF_mp(grfp, fasta_files_list, t, length)
        collect_results(genome_name, fasta_files_list)
        subprocess.Popen(["rm", "-rf", "GRFmite_mp"])
