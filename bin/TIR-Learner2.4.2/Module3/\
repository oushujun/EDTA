import pandas as pd
import numpy as np
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
parser.add_argument("-s", "--species", help="One of the following: Maize , Rice or others", required=True)
args = parser.parse_args()#pylint: disable=invalid-name

genome_file=args.genomeFile
genome_Name = args.genomeName
path=args.path
t=args.processer
dir=args.currentD


targetDir=dir+"/"+genome_Name+"/"

os.chdir(targetDir)

spliter="-+-"

species=args.species

if species == "others":
    species="Maize"

def en5(file):
    f=open(file,"r+")
    lines=f.readlines()
    if len(lines)>1:
        topre = pd.read_csv(file, header=0)
        toprex = topre.drop(['ID'], axis=1)
        filename = path+"/Module3/"+species+'_model.sav'
        ensemblecf = joblib.load(filename)
        pre = ensemblecf.predict(toprex)
        p = pd.DataFrame(pre,columns=['prediction'])
        r = pd.DataFrame(topre['ID']).join(p)
        r.to_csv(file+spliter+"predi.csv",index=None)
    else:
        print("No candidate found in %s, skip this"%(file))    



if __name__ == '__main__':
    files=os.listdir(".")
    files=[i for i in files if i.split(spliter)[-1]=="toPre.csv"]
    pool = multiprocessing.Pool(int(t))
    pool.map(en5,files)
    pool.close()
    pool.join()
    

