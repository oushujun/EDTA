import os
import subprocess
import multiprocessing as mp

import prog_const


def blast_ref_lib_in_genome_file(genome_db, genome_name, ref_lib, ref_lib_file_path, cpu_cores):
    out = genome_name + prog_const.spliter + "blast" + prog_const.spliter + ref_lib
    blast = (f"blastn -max_hsps 5 -perc_identity 80 -qcov_hsp_perc 100 -query \"{ref_lib_file_path}\" "
             f"-db \"{genome_db}\" -num_threads {int(cpu_cores)} -outfmt '6 qacc sacc length pident gaps mismatch "
             f"qstart qend sstart send evalue qcovhsp' -out \"{out}\" 2>/dev/null")
    # Find where is the RefLib in the genome database
    subprocess.Popen(blast, shell=True).wait()


def blast_de_novo_result_in_ref_lib(file_name, ref_lib, ref_lib_file_path, cpu_cores):
    out = file_name + prog_const.spliter + "blast" + prog_const.spliter + ref_lib
    blast = (f"blastn -max_hsps 5 -perc_identity 80 -qcov_hsp_perc 80 -query \"{file_name}\" "
             f"-subject \"{ref_lib_file_path}\" -num_threads {int(cpu_cores)} -outfmt '6 qseqid sseqid length pident "
             f"gaps mismatch qstart qend sstart send evalue qcovhsp' -out \"{out}\" 2>/dev/null")
    # Find whether there is any predicted TIR (GRFmite and TIRvish) inside the refLib
    subprocess.Popen(blast, shell=True).wait()

# def collect_blast_result():
#     df = pd.concat([pd.read_csv(f, header=None) for f in
#                     [f for f in os.listdir(".") if f.endswith("_RefLib") and os.path.getsize(f) != 0]])
#     csv_QUOTE_NONE = 3
#     df.to_csv("tem_blastResult", header=False, index=False, quoting=csv_QUOTE_NONE)


def blast_genome_file(TIRLearner_instance):
    print("Module 1, Step 1: Blast reference library in genome file")
    genome_db = TIRLearner_instance.genome_file_path + prog_const.spliter + "db"
    mkDB = (f"makeblastdb -in {TIRLearner_instance.genome_file_path} -out {genome_db} "
            f"-parse_seqids -dbtype nucl 2>/dev/null")
    subprocess.Popen(mkDB, shell=True).wait()

    mp_args_list = [(genome_db, TIRLearner_instance.genome_name,
                     ref_lib, os.path.join(prog_const.ref_lib_dir_path, ref_lib),
                     TIRLearner_instance.cpu_cores)
                    for ref_lib in prog_const.ref_lib_file_dict[TIRLearner_instance.species]]

    print(mp_args_list)

    with mp.Pool(int(TIRLearner_instance.cpu_cores)) as pool:
        pool.starmap(blast_ref_lib_in_genome_file, mp_args_list)

    subprocess.Popen(["find", ".", "-name", f"*{prog_const.spliter}db*", "-delete"])  # remove blast db files


def blast_de_novo_result(TIRLearner_instance):
    print("Module 2, Step 6: Blast GRF and TIRvish result in reference library")
    mp_args_list = [(TIRLearner_instance.processed_de_novo_result_file_name,
                     ref_lib, os.path.join(prog_const.ref_lib_dir_path, ref_lib),
                     TIRLearner_instance.cpu_cores)
                    for ref_lib in prog_const.ref_lib_file_dict[TIRLearner_instance.species]]

    with mp.Pool(int(TIRLearner_instance.cpu_cores)) as pool:
        pool.starmap(blast_de_novo_result_in_ref_lib, mp_args_list)
