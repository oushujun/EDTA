import pandas as pd
import numpy as np
#from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier, VotingClassifier
#from sklearn.neighbors import KNeighborsClassifier
#from sklearn.model_selection import cross_val_score
#from sklearn.neural_network import MLPClassifier
#from sklearn import preprocessing as pp
from sklearn.externals import joblib

import warnings
import time
import os
from sklearn.tree import DecisionTreeClassifier
import argparse
import multiprocessing
from multiprocessing import Pool
from Bio import SeqIO
import subprocess
# remove warnings
def warn(*args, **kwargs):
    pass

warnings.warn = warn
warnings.filterwarnings("ignore")


seed=7

parser = argparse.ArgumentParser()#pylint: disable=invalid-name
parser.add_argument("-g", "--genomeFile", help="Genome file in fasta format", required=True)

parser.add_argument("-name", "--genomeName", help="Genome Name", required=True)
parser.add_argument("-p", "--path", help="Source code path", required=True)
parser.add_argument("-t", "--processer", help="Number of processer", required=True)
parser.add_argument("-d", "--currentD", help="Path of current directory", required=True)
args = parser.parse_args()#pylint: disable=invalid-name

genome_file=args.genomeFile
genome_Name = args.genomeName
path=args.path
t=args.processer
dir=args.currentD


targetDir=dir+"/"+genome_Name+"/"

spliter="-+-"


def en5(file):
    topre = pd.read_csv(file, header=0)
    toprex = topre.drop(['ID'], axis=1)
    filename = path+"/Module3/"+'Rice_model.sav'
    ensemblecf = joblib.load(filename)
    pre = ensemblecf.predict(toprex)
    p = pd.DataFrame(pre,columns=['prediction'])
    r = pd.DataFrame(topre['ID']).join(p)
    r.to_csv(file+spliter+"predi.csv",index=None)


#arglist(genomefile,selectfile)
def GetFastaFromFile(argList):
    genomedb=argList[0]
    prediction=argList[1]
    f = pd.read_csv(prediction, header=0, sep=",") #shujun
  #  f = pd.read_table(prediction, header=0, sep=",")
    for i in ["DTA","DTC","DTH","DTM","DTT"]:
        sub=f.loc[f["prediction"]==i]
        preNamelist=list(sub["ID"])
        outname=prediction[:-4]+spliter+i+".fa"
        o = open(outname,"w")
        for line in preNamelist:
            infor=line
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
            o.write(">"+line+spliter+entry+spliter+str(p1)+spliter+str(p2)+spliter+str(p_start)+spliter+str(p_end) + "\n" + str(out_seq1) + "\n")
        o.close()

# #





if __name__ == '__main__':
    os.chdir(targetDir)
    files=os.listdir(".")
    files=[i for i in files if i.split(spliter)[-1]=="toPre.csv"]
    pool = multiprocessing.Pool(int(t))
    pool.map(en5,files)
    pool.close()
    pool.join()
    
    files=os.listdir(".")
    fastafile=[i for i in files if i[-10:]=="nonHomo.fa"]
    lists=[[genome_file+spliter+"db", j[:-3]+spliter+"toPre.csv"+spliter+"predi.csv"] for j in fastafile]
    pool = multiprocessing.Pool(int(t))
    pool.map(GetFastaFromFile,lists)
    pool.close()
    pool.join()
    



