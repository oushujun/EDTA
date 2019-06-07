from Bio import SeqIO
import multiprocessing
from multiprocessing import Pool
import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser()#pylint: disable=invalid-name
parser.add_argument("-g", "--genomeFile", help="Genome file in fasta format", required=True)
parser.add_argument("-name", "--genomeName", help="Genome Name", required=True)
parser.add_argument("-p", "--path", help="Source code path", required=True)
parser.add_argument("-t", "--processer", help="Number of processer", required=True)
parser.add_argument("-d", "--currentD", help="Path of current directory", required=True)
parser.add_argument("-grfp", "--GRFpath", help="Path of GRF program", required=True)
args = parser.parse_args()#pylint: disable=invalid-name

genome_file = args.genomeFile
genome_Name = args.genomeName
path=args.path
t=args.processer
dir=args.currentD
grfP=args.GRFpath
spliter="-+-"
targetDir=dir+"/"+genome_Name+"/"



grf="%s/grf-main -i %s -o %s -c 1 -t %s -p 20 --min_space 10 " \
        "--max_space 10000 --max_indel 0 --min_tr 10 --min_spacer_len 10 --max_spacer_len 10000"%(grfP,genome_file,genome_Name,int(t))
os.system(grf)



