import os
import subprocess
import multiprocessing as mp

spliter = "-+-"
TIR_types = ("DTA", "DTC", "DTH", "DTM", "DTT")


def blast_reference(file_name, refLib, path):
    subject = os.path.join(path, "RefLib", refLib)
    # subject = path + "/RefLib/" + refLib
    out = file_name + spliter + "blast" + spliter + refLib
    blast = (f"blastn -max_hsps 5 -perc_identity 80 -qcov_hsp_perc 80 -query \"{file_name}\" -subject \"{subject}\" "
             "-outfmt '6 qseqid sseqid length pident gaps mismatch qstart qend sstart send evalue qcovhsp' "
             f"-out \"{out}\" 2>/dev/null")
    # Find whether there is any predicted TIR (GRFmite) inside the refLib
    subprocess.Popen(blast, shell=True).wait()


def execute(args):
    print("Module 2, Step 3A: GRF result blast reference sequences")
    genome_name = args[1]
    processedGRFmite_file_name = f"{genome_name}{spliter}processedGRFmite.fa"
    path = args[2]
    t = args[3]
    species = args[4]

    ref_list = [f"{species}_{TIR_type}_RefLib" for TIR_type in TIR_types]
    mp_args_list = [(processedGRFmite_file_name, ref, path) for ref in ref_list]

    with mp.Pool(int(t)) as pool:
        pool.starmap(blast_reference, mp_args_list)
