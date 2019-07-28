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
parser.add_argument("-l", "--TIR_length", help="Max Length of TIR", required=True, default=10000)
args = parser.parse_args()#pylint: disable=invalid-name

genome_file = args.genomeFile
genome_Name = args.genomeName
path=args.path
t=args.processer
dir=args.currentD
spliter="-+-"
targetDir=dir+"/"+genome_Name+"/"
length=args.TIR_length
os.chdir(targetDir)


print("GT starts indexing")
pathcode=path+"/Module2/"
Index=pathcode+"gt-1.5.10-Linux_x86_64-64bit-complete/bin/gt suffixerator -db %s -indexname %s -tis -suf -lcp -des -ssp -sds -dna -mirrored"%(genome_file,genome_Name)
print(Index)
os.system(Index)
print("Indexing Finished")


print("GT starts searching for TIR candidates")

RunGT=pathcode+"gt-1.5.10-Linux_x86_64-64bit-complete/bin/gt tirvish -index %s -seed 20 -mintirlen 10 -maxtirlen 1000 -mintirdist 10 -maxtirdist %s -similar 80 -mintsd 2 -maxtsd 11 -vic 13 -seqids 'yes' > %s_TIRvish.gff"%(genome_Name,length,genome_Name)
os.system(RunGT)

print("GT finished")

