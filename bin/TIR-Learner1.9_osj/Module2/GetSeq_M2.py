import subprocess
import multiprocessing
from multiprocessing import Pool
import argparse
import os
from Bio import SeqIO
import pandas as pd
import glob #shujun


parser = argparse.ArgumentParser()#pylint: disable=invalid-name
parser.add_argument("-g", "--genomeFile", help="Genome file in fasta format", required=True)
parser.add_argument("-name", "--genomeName", help="Genome Name", required=True)
parser.add_argument("-p", "--path", help="Source code path", required=True)
parser.add_argument("-d", "--currentD", help="Path of current directory", required=True)
parser.add_argument("-t", "--processer", help="Number of processer", required=True)

args = parser.parse_args()#pylint: disable=invalid-name

genome_file = args.genomeFile
genome_Name = args.genomeName
path=args.path
t=args.processer
dir=args.currentD

targetDir=dir+"/"+genome_Name+"/"
spliter="-+-"


#arglist(genomefile,selectfile)
def GetFastaFromFile(argList):
    genomedb=argList[0]
    selectfile=argList[1]
    f=open(selectfile,"r+")
    outname=selectfile+".fa"
    o = open(outname,"w")
    lines=f.readlines()
    if len(lines)!=0:
        for line in lines:
            infor=line.split("\t")[0]
            entry=infor.split(":")[0]
            p1 = int(infor.split(":")[1])
            p2 = int(infor.split(":")[2])
            if (p1 < p2):
                if (p1>200):
                    p_start = p1-200
                else:
                    p_start=1
                p_end = p2 + 200

            else:
                if (p2 > 200):
                    p_start = p2-200
                else:
                    p_start = 1
                p_end = p1+200
            seq1 = subprocess.check_output("blastdbcmd -db '%s' -entry '%s' -range '%s'-'%s'" % (genomedb, entry, int(p_start), int(p_end)), shell=True)
            seq1 = seq1.decode("utf-8")
            out_seq1 = ''
            split = seq1.split('\n')
            for sp in split:
                if any([i.isdigit() for i in sp]):
                    continue
                out_seq1 += sp
            o.write(">"+line.split("\t")[0]+spliter+entry+spliter+str(p1)+spliter+str(p2)+spliter+str(p_start)+spliter+str(p_end) + "\n" + str(out_seq1) + "\n")
        f.close()
        o.close()

# #
# #
def getNonHomo(argList):
    allSeq = argList[0]
    homoseq=argList[1]
    records_all=list(SeqIO.parse(allSeq,"fasta"))
    homoseq=list(SeqIO.parse(homoseq,"fasta"))
    homoID=[rec.id.split("_")[0] for rec in homoseq]
    SeqIO.write((rec for rec in records_all if rec.id not in homoID), allSeq[:-4]+"_nonHomo.fa","fasta")


if __name__ == '__main__':
    os.chdir(targetDir)
  #  fileList=os.listdir(".")
  #  fileList=[i for i in fileList if i[-2:]=="80"]
    fileList=glob.glob('./*GRFmite.fa*80') #shujun

    genomedb=genome_file+spliter+"db"
    argList=[[genomedb,selectfile] for selectfile in fileList]
    pool = multiprocessing.Pool(int(t))
    pool.map(GetFastaFromFile,argList)
    pool.close()
    pool.join()

    names=pd.read_csv("%sContig.name"%(genome_Name+spliter),header=None)
  #  names=pd.read_table("%sContig.name"%(genome_Name+spliter),header=None)
    l_name=list(names[0])
    for name in l_name:
        cat= "cat *%s*80.fa > %s"%(name,genome_Name+spliter+name+spliter+"homo.fa")
        os.system(cat)
    l_list=[["%sGRFmite.fa-+-p"%(genome_Name+spliter+name+spliter), "%s"%(genome_Name+spliter+name+spliter+"homo.fa")] for name in l_name]
    pool = multiprocessing.Pool(int(t))
    pool.map(getNonHomo,l_list)
    pool.close()
    pool.join()
  
