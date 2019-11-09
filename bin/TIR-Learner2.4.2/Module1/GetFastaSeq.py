import multiprocessing
from multiprocessing import Pool
import argparse
import os
from Bio import SeqIO
import numpy as np
import subprocess

parser = argparse.ArgumentParser()#pylint: disable=invalid-name
parser.add_argument("-name", "--genomeName", help="Genome Name", required=True)
parser.add_argument("-p", "--path", help="Source code path", required=True)
parser.add_argument("-t", "--processer", help="Number of processer", required=True)
parser.add_argument("-d", "--currentD", help="Path of current directory", required=True)
parser.add_argument("-g", "--genomeFile",help="Path to the genome file", required=True)
args = parser.parse_args()#pylint: disable=invalid-name

genome_Name = args.genomeName
path=args.path
t=args.processer
dir=args.currentD
genomeFile=args.genomeFile

targetDir=dir+"/"+genome_Name+"/"
os.chdir(targetDir)
spliter="-+-"


def GetListFromFile(file):
    listName=file+spliter+"200.list" #shujun
    records=list(SeqIO.parse(file,"fasta"))
    f=open(file,"r+")
    lines=f.readlines()
    o = open(listName,"a+")
    for line in lines: #shujun
        p1 = int(line.split("\t")[8])
        p2 = int(line.split("\t")[9])
        entry=line.split("\t")[1]
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
        o.write(genome_Name+spliter+entry+spliter+str(p1)+spliter+str(p2)+spliter+str(p_start)+spliter+str(p_end)+ "\t" + entry + ":" + str(p_start) + ".." + str(p_end) + "\n")
    o.close() #shujun


def GetFastaFromList(argList): #shujun
    genomeFile=argList[0]
    listName=argList[1]
    outName=argList[2]
    get_seq = "perl %s/Module3_New/call_seq_by_list2.pl %s -C %s -header 1 -out %s" % (path, listName, genomeFile, outName)
    subprocess.run(['/bin/bash', '-c', get_seq])
    clean_head = "perl -i -nle 's/^>.*\|/>/; print $_' %s" % (outName)
    subprocess.run(['/bin/bash', '-c', clean_head])


if __name__ == '__main__':
    files=os.listdir(".")
    fileList=[i for i in files if i[-3:]=="csv"]
    argList=[[genomeFile, fileList[i]+spliter+"200.list",fileList[i][:-4]+".fa"] for i in range(0,len(fileList))] #shujun
    pool = multiprocessing.Pool(int(t))
    pool.map(GetListFromFile,fileList) #shujun
    pool.map(GetFastaFromList,argList) #shujun
    pool.close()
    pool.join()
