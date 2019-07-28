import pandas as pd
from Bio import SeqIO
from itertools import product
import regex as re
import multiprocessing
from multiprocessing import Pool
import argparse
import os

def getK_mers(k):
    N="ATCG"
    L=product(N,repeat=k)
    t=[i for i in list(L)]
    k_mers=[]
    for i in t:
        s=''
        for j in i:
            s+=j
        k_mers.append(s)
    return k_mers


def getFeatureList(n_kmer):
    featureList=["ID"]
    for i in range(1,n_kmer):
        l=getK_mers(i)
        featureList+=l
    return featureList


def getTrainingset(rec):
    featureList=getFeatureList(5)
    Dic=[]
    Dic.append(str(rec.id))
    for i in featureList[1:]:
        seq=str(rec.seq)[20:-20].upper()
        Dic.append(len(re.findall(i, seq, overlapped=True)))
    rec_df=pd.DataFrame([Dic],columns=featureList)
    v=""
    for k in Dic[:-1]:
        v+=str(k)+","
    v+=str(Dic[-1])
    return Dic


parser = argparse.ArgumentParser()#pylint: disable=invalid-name
parser.add_argument("-name", "--genomeName", help="Genome Name", required=True)
parser.add_argument("-p", "--path", help="Source code path", required=True)
parser.add_argument("-t", "--processer", help="Number of processer", required=True)
parser.add_argument("-d", "--currentD", help="Path of current directory", required=True)
args = parser.parse_args()#pylint: disable=invalid-name

genome_Name = args.genomeName
path=args.path
t=args.processer
dir=args.currentD


targetDir=dir+"/"+genome_Name+"/"
os.chdir(targetDir)

spliter="-+-"


if __name__ == '__main__':
    files=os.listdir(".")
    files=[i for i in files if i[-10:]=="nonHomo.fa"]
    for file in files:
        out = open("%stoPre.csv"%(file[:-3]+spliter), "w")
       # out = open("%stoPre.csv"%(file[:-3]+spliter), "a+")
        s = ""
        featureList = getFeatureList(5)
        for j in featureList[:-1]:
            s += j + ","
        s += featureList[-1]
        out.write(s + "\n")
        records=list(SeqIO.parse(file,"fasta"))
        
        pool = multiprocessing.Pool(int(t))
        d = pool.map(getTrainingset, records)
        pool.close()
        pool.join()
        for single in d:
            data = ""
            for init in single[:-1]:
                data += str(init) + ","
            data += str(single[-1])
            out.write(data + "\n")
        out.close()


