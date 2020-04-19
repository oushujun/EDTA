import multiprocessing
from multiprocessing import Pool
import argparse
import os
from Bio import SeqIO
import numpy as np
import subprocess
import glob
import pandas as pd
from os import path
import os.path

parser = argparse.ArgumentParser()#pylint: disable=invalid-name
parser.add_argument("-name", "--genomeName", help="Genome Name", required=True)
parser.add_argument("-p", "--path", help="Source code path", required=True)
parser.add_argument("-t", "--processer", help="Number of processer", required=True)
parser.add_argument("-d", "--currentD", help="Path of current directory", required=True)
parser.add_argument("-g", "--genomeFile",help="Path to the genome file", required=True)
args = parser.parse_args()#pylint: disable=invalid-name

genome_Name = args.genomeName
pathi=args.path
t=args.processer
dir=args.currentD
genomeFile=args.genomeFile

targetDir=dir+"/"+genome_Name+"/"
os.chdir(targetDir)
spliter="-+-"


def GetListFromFile(file):
    f=open(file,"r+")
    lines=f.readlines()
    listName=file+spliter+"200.list" #shujun
    o = open(listName,"a+")
    if len(lines)!=0:
        for line in lines: #shujun
            infor=line.split("\t")[0]
            entry=infor.split(":")[0]
            p1 = int(infor.split(":")[1])
            p2 = int(infor.split(":")[2])
            fam=infor.split("_")[-1]
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
            o.write(line.split("\t")[0]+spliter+entry+spliter+str(p1)+spliter+str(p2)+spliter+str(p_start)+spliter+str(p_end)+ "\t" + entry + ":" + str(p_start) + ".." + str(p_end) + "\n")
        o.close() #shujun


def GetFastaFromList(argList): #shujun
    genomeFile=argList[0]
    listName=argList[1]
    outName=argList[2]
    get_seq = "perl %s/Module3_New/call_seq_by_list2.pl %s -C %s -header 1 -out %s" % (pathi, listName, genomeFile, outName)
    subprocess.run(['/bin/bash', '-c', get_seq])
    clean_head = "perl -i -nle 's/^>.*\|/>/; print $_' %s" % (outName)
    subprocess.run(['/bin/bash', '-c', clean_head])


def getNonHomo(argList):
    allSeq = argList[0]
    homoseq=argList[1]
    if (path.isfile(homoseq)==True and path.isfile(allSeq)==True):
        records_all=list(SeqIO.parse(allSeq,"fasta"))
        homoseq=list(SeqIO.parse(homoseq,"fasta"))
        homoID=[rec.id.split(spliter)[0] for rec in homoseq]
        SeqIO.write((rec for rec in records_all if rec.id not in homoID), allSeq[:-4]+"_nonHomo.fa","fasta")


if __name__ == '__main__':
    files=os.listdir(".")
    fileList=glob.glob('./*GRFmite.fa*80')
    argList=[[genomeFile,fileList[i]+spliter+"200.list",fileList[i]+".fa"] for i in range(0,len(fileList))] #shujun
    pool = multiprocessing.Pool(int(t))
    pool.map(GetListFromFile,fileList) #shujun
    pool.map(GetFastaFromList,argList) #shujun
    pool.close()
    pool.join()
   
    names=pd.read_csv("%sContig.name"%(genome_Name+spliter),header=None)
    l_name=list(names[0].astype(str))
    for name in l_name:
       # cat= "cat *%s*80.fa > %s"%(name,genome_Name+spliter+name+spliter+"homo.fa")
        cat="for i in *%s*80.fa; do cat $i; done > %s"%(name,genome_Name+spliter+name+spliter+"homo.fa") #shujun
        os.system(cat)
    l_list=[["%sGRFmite.fa-+-p"%(genome_Name+spliter+name+spliter), "%s"%(genome_Name+spliter+name+spliter+"homo.fa")] for name in l_name]
    pool = multiprocessing.Pool(int(t))
    pool.map(getNonHomo,l_list)
    pool.close()
    pool.join()
