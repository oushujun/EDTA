import os
import subprocess
import multiprocessing as mp

import prog_const
spliter = prog_const.spliter
TIR_types = prog_const.TIR_types
path = prog_const.program_root_dir



def blast_reference(genome_db, genome_name, ref_lib, path, t):
    query = os.path.join(path, "RefLib", ref_lib)
    print(query)
    # query = path + "/RefLib/" + refLib
    out = genome_name + spliter + "blast" + spliter + ref_lib
    blast = f"blastn -max_hsps 5 -perc_identity 80 -qcov_hsp_perc 100 -query \"{query}\" -db \"{genome_db}\" " + \
            f"-num_threads {int(t)} -outfmt '6 qacc sacc length pident gaps mismatch qstart qend sstart " + \
            f"send evalue qcovhsp' -out \"{out}\" 2>/dev/null"
    # Find where is the RefLib in the genome database
    subprocess.Popen(blast, shell=True).wait()


# def collect_blast_result():
#     df = pd.concat([pd.read_csv(f, header=None) for f in
#                     [f for f in os.listdir(".") if f.endswith("_RefLib") and os.path.getsize(f) != 0]])
#     csv_QUOTE_NONE = 3
#     df.to_csv("tem_blastResult", header=False, index=False, quoting=csv_QUOTE_NONE)


def execute(TIRLearner_instance):
    print("Module 1, Step 1: Blast Genome against Reference Library")
    genome_file = TIRLearner_instance.genome_file
    genome_name = TIRLearner_instance.genome_name
    species = TIRLearner_instance.species
    t = TIRLearner_instance.cpu_cores

    genome_db = genome_file + spliter + "db"
    mkDB = f"makeblastdb -in {genome_file} -out {genome_db} -parse_seqids -dbtype nucl 2>/dev/null"
    subprocess.Popen(mkDB, shell=True).wait()

    ref_lib_list = [f"{species}_{TIR_type}_RefLib" for TIR_type in TIR_types]
    mp_args_list = [(genome_db, genome_name, ref_lib, path, t) for ref_lib in ref_lib_list]

    with mp.Pool(int(t)) as pool:
        pool.starmap(blast_reference, mp_args_list)

    # remove blast db files
    # subprocess.Popen(["rm", "-f", f"../*{spliter}db*"])
    # subprocess.Popen(["find", "../", "-name", f"*{spliter}db*", "-delete"])
    subprocess.Popen(["find", ".", "-name", f"*{spliter}db*", "-delete"])
