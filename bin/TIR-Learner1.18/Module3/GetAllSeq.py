import subprocess
import multiprocessing
from multiprocessing import Pool
import argparse
import os
from Bio import SeqIO
import pandas as pd

parser = argparse.ArgumentParser()#pylint: disable=invalid-name
parser.add_argument("-g", "--genomeFile", help="Genome file in fasta format", required=True)
parser.add_argument("-name", "--genomeName", help="Genome Name", required=True)
parser.add_argument("-p", "--path", help="Source code path", required=True)
parser.add_argument("-t", "--processer", help="Number of processer", required=True)
parser.add_argument("-d", "--currentD", help="Path of current directory", required=True)
args = parser.parse_args()#pylint: disable=invalid-name

genome_file = args.genomeFile
genome_Name = args.genomeName
path=args.path
t=args.processer
dir=args.currentD

targetDir=dir+"/"+genome_Name+"/"

spliter="-+-"


os.chdir(targetDir)

#arglist(genomefile,selectfile)
def GetFastaFromFile(argList):
    genomedb=argList[0]
    line=argList[1]
    outname=argList[2]
    gName=argList[3]
   # o = open(outname,"w")
    o = open(outname,"a+") #shujun
    entry=line.split("\t")[0]
    p1 = int(line.split("\t")[3])
    p2 = int(line.split("\t")[4])
    family=line.split("\t")[2]
    infor=line.split("\t")[8]
    seq1 = subprocess.check_output("blastdbcmd -db '%s' -entry '%s' -range '%s'-'%s'" % (genomedb, entry, int(p1), int(p2)), shell=True)
    seq1 = seq1.decode("utf-8")
    out_seq1 = ''
    split = seq1.split('\n')
    for sp in split:
        if any([i.isdigit() for i in sp]):
            continue
        out_seq1 += sp
   # o.write(">%s_%s_%s_%s_%s_%s" % (gName,entry, p1, p2,family,infor) + "\n" + str(out_seq1) + "\n")
    o.write(">%s_%s_%s_%s_%s_%s" % (gName,entry, p1, p2,family,infor) + str(out_seq1) + "\n") #shujun


if __name__ == '__main__':
    genomedb=genome_file+spliter+"db"
#    f="%s_FinalAnn.gff3"%(genome_Name+"_combine"+"/"+genome_Name)
    f="%s_FinalAnn.gff3"%(genome_Name)
    outName=f[0:-5]+".fa"
    file=open(f,"r+")
    lines=file.readlines()
    argList=[[genomedb,i,outName,genome_Name] for i in lines]
    pool = multiprocessing.Pool(int(t))
    pool.map(GetFastaFromFile,argList)
    pool.close()
    pool.join()

#    f="%s_FinalAnn_Clint.gff3"%(genome_Name)
#    outName=f[0:-5]+".fa"
#    file=open(f,"r+")
#    lines=file.readlines()
#    argList=[[genomedb,i,outName,genome_Name] for i in lines]
#    pool = multiprocessing.Pool(16)
#    pool.map(GetFastaFromFile,argList)
#    pool.close()
#    pool.join()



