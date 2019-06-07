from Bio import SeqIO
import multiprocessing
from multiprocessing import Pool
import pandas as pd
import argparse
import regex as re
import os

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
spliter="-+-"
targetDir=dir+"/"+genome_Name+"/"



def getContigNames(genomeFile,genomeName):
    records=list(SeqIO.parse(genomeFile,"fasta"))
    l=[rec.id for rec in records]
    f=pd.DataFrame(l)
    f.to_csv("%sContig.name"%(genomeName+spliter),header=None,index=None)



def SplitFasta(contigFile,FastaFile,genomeName):
    records=list(SeqIO.parse(FastaFile,"fasta"))
    names=pd.read_table(contigFile,header=None)
    l_name=list(names[0])
    for name in l_name:
        SeqIO.write((rec for rec in records if rec.id.split(":")[0]==name),"%sGRFmite.fa"%(genomeName+spliter+name+spliter),"fasta")



def TArepeats(s):
    t=s.upper().count("T")
    a=s.upper().count("A")
    ta=t+a
    if (ta>len(s)*0.7):
        return True
    else:
        return False


def checkN(s):
    n=s.upper().count("N")
    if n>0:
        return True
    else:
        return False

def checkNPer(s):
    n=s.upper().count("N")
    p=n/len(s)
    if p>=0.20:
        return True
    else:
        return False

def findDigitsSum(string):
    pattern = '(\d+)'
    l = re.findall(pattern,string)
    return sum([int(i) for i in l])


def getSeqID(file):
    remove=[]
    records=list(SeqIO.parse(file,"fasta"))
    for rec in records:
        s=str(rec.seq)
        tirLen = findDigitsSum(rec.id.split(":")[-2])
        tir=str(rec.seq)[0:tirLen]
        print(rec.id+" "+str(tirLen))
        if (TArepeats(s)==True or checkNPer(s)==True or TArepeats(tir)==True or checkN(tir)==True or len(s)<50):
            remove.append(rec.id)
    records=list(SeqIO.parse(file,"fasta"))
    l=[]
    for rec in records:
        tsd=rec.id.split(":")[-1]
        if (len(tsd)>6 or tsd=="TAA" or tsd=="TTA" or tsd=="TA" or str(rec.seq)[0:4]=="CACT") or str(rec.seq)[0:4]=="GTGA":
            l.append(rec.id)
    SeqIO.write((rec for rec in records if rec.id in l and rec.id not in remove), file+"_p","fasta")
    


if __name__ == '__main__':
    os.chdir("./%s"%(genome_Name))
    getContigNames(genome_file, genome_Name)
    SplitFasta("%s_Contig.name"%(genome_Name), "candidate.fasta", genome_Name)
    Allfiles=os.listdir(".")
    files=[i for i in Allfiles if i[-10:]=="GRFmite.fa"]
    pool = multiprocessing.Pool(16)
    pool.map(getSeqID,files)
    pool.close()
    pool.join()





