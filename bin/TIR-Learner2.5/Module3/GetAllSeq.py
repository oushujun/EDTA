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

genomeFile = args.genomeFile
genome_Name = args.genomeName
path=args.path
t=args.processer
dir=args.currentD

targetDir=dir+"/"+genome_Name+"/"

spliter="-+-"


os.chdir(targetDir)


def GetListFromFile(file):
    f=open(file,"r+")
    lines=f.readlines()
    listName=file+spliter+".list" #shujun
    o = open(listName,"a+")
    if len(lines)!=0:
        for line in lines: #shujun
            entry=line.split("\t")[0]
            p1 = int(line.split("\t")[3])
            p2 = int(line.split("\t")[4])
            family=line.split("\t")[2]
            infor=line.split("\t")[8]
            o.write("%s_%s_%s_%s_%s_%s" % (genome_Name,entry, p1, p2,family,infor[0:-1])+ "\t" + entry + ":" + str(p1) + ".." + str(p2) + "\n")
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
    f="%s_FinalAnn.gff3"%(genome_Name)
    fileList=["%s_FinalAnn.gff3"%(genome_Name),"%s_FinalAnn_filter.gff3"%(genome_Name)]
    argList=[[genomeFile,fileList[i]+spliter+".list",fileList[i][0:-5]+".fa"] for i in range(0,len(fileList))] #shujun
    pool = multiprocessing.Pool(int(t))
    pool.map(GetListFromFile,fileList) #shujun
    pool.map(GetFastaFromList,argList) #shujun
    pool.close()
    pool.join()


