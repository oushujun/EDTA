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
parser.add_argument("-l", "--TIR_length", help="Max Length of TIR", required=True, default=5000)
args = parser.parse_args()#pylint: disable=invalid-name

genome_file = args.genomeFile
genome_Name = args.genomeName
path=args.path
t=args.processer
dir=args.currentD
grfP=args.GRFpath
spliter="-+-"
targetDir=dir+"/"+genome_Name+"/"
length=args.TIR_length
tempDir=dir+"/"+"temp"+"/"


def getContigNames(genomeFile,genomeName):
    records=list(SeqIO.parse(genomeFile,"fasta"))
    l=[rec.id for rec in records]
    f=pd.DataFrame(l)
    f.to_csv(targetDir+"%sContig.name"%(genomeName+spliter),header=None,index=None)

def SplitFasta(contigFile,genomeFile,genomeName):
    records=list(SeqIO.parse(genomeFile,"fasta"))
    names=pd.read_csv(contigFile,header=None) #shujun
   # names=pd.read_table(contigFile,header=None)
    l_name=list(names[0].astype(str))
    for name in l_name:
        SeqIO.write((rec for rec in records if rec.id.split(":")[0]==name),"%s.fasta"%(targetDir+genomeName+spliter+name),"fasta")

os.chdir(tempDir)
files=os.listdir(".")
p_files=[i for i in files if i.split(spliter)[0]==genome_Name and i.split(spliter)[-1]=="p" and i.split(spliter)[-2]=="GRFmite.fa"]
grf_files=[i for i in files if i.split(spliter)[0]==genome_Name and i[-10:]=="GRFmite.fa"]

getContigNames(genome_file, genome_Name)

if (len(p_files)>0):
    n=0
    for i in p_files:
        f=open(i,"r+")
        lines=f.readlines()
        n+=len(lines)
    if (n>0):
         print("Copy processed GRF output from temp to %s"%(targetDir))
         cp= "cp %s*GRFmite.fa-+-p %s"%(genome_Name+spliter,targetDir)
         os.system(cp)
         rm="rm  %s*GRFmite.fa > /dev/null 2>&1"%(genome_Name+spliter) #shujun
         print("Processed GRF output files exist, remove raw GRF output files")
         os.system(rm)
elif (len(grf_files)>0):
    print("Copy raw GRF output from temp to %s"%(targetDir))
    cp= "cp %s*GRFmite.fa %s"%(genome_Name+spliter,targetDir)
    os.system(cp)
else:
    print("Running GRF")
    getContigNames(genome_file, genome_Name)
    os.chdir(targetDir)
    SplitFasta("%sContig.name"%(genome_Name+spliter),genome_file,genome_Name)
    
    files=os.listdir(".")
    files=[i for i in files if i.split(spliter)[0]==genome_Name and i[-6:]==".fasta"]
    
    for file in files:
        records=list(SeqIO.parse(file,"fasta"))
        if (len(str(records[0].seq))>int(length)+500):
            grf="%s/grf-main -i %s -o %s -c 1 -t %s -p 20 --min_space 10 " \
            "--max_space %s --max_indel 0 --min_tr 10 --min_spacer_len 10 --max_spacer_len %s"%(grfP,file,file.split(spliter)[1][:-6],int(t),int(length),int(length))
            os.system(grf)
    
    

