import multiprocessing
from multiprocessing import Pool
import argparse
import os
from Bio import SeqIO
import numpy as np
import subprocess
import pandas as pd

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


def GetListFromFile(arglist):
    lines=arglist[0]
    listName=arglist[1] #shujun
    o = open(listName,"a+")
    if len(lines)!=0:
        for line in lines: #shujun
            infor=line
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
            o.write(line+spliter+entry+spliter+str(p1)+spliter+str(p2)+spliter+str(p_start)+spliter+str(p_end)+ "\t" + entry + ":" + str(p_start) + ".." + str(p_end) + "\n")
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
    fastafile=[i for i in files if i[-10:]=="nonHomo.fa"]
    prefiles=[j[:-3]+spliter+"toPre.csv"+spliter+"predi.csv" for j in fastafile]
    for prediction in prefiles:
        if os.path.exists(prediction)==True:
            f = pd.read_csv(prediction, header=0, sep=",")
            for i in ["DTA","DTC","DTH","DTM","DTT"]:
                sub=f.loc[f["prediction"]==i]
                preNamelist=list(sub["ID"])
                argList1=[[preNamelist,prediction[:-4]+spliter+i+".list"]]
                argList=[[genomeFile,prediction[:-4]+spliter+i+".list",prediction[:-4]+spliter+i+".fa"]] #shujun
                pool = multiprocessing.Pool(int(t))
                pool.map(GetListFromFile,argList1) #shujun
                pool.map(GetFastaFromList,argList) #shujun
                pool.close()
                pool.join()
