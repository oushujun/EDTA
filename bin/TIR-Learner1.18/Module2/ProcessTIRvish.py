from Bio import SeqIO
import multiprocessing
from multiprocessing import Pool
import pandas as pd
import argparse
import regex as re
import os
import subprocess

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

os.chdir(targetDir)


def getGff(file,out):
    f=open(file,"r")
    o=open(out,"w")
    lines=f.readlines()
    i=0
    while (i<len(lines)):
        if (lines[i]=="###\n"):
            cluster=lines[i-6:i]
            tsd1=cluster[1].split("\t")[3]+"_"+cluster[1].split("\t")[4]
            tir1=cluster[3].split("\t")[3]+"_"+cluster[3].split("\t")[4]
            tir2=cluster[4].split("\t")[3]+"_"+cluster[4].split("\t")[4]
            tsd2=cluster[5].split("\t")[3]+"_"+cluster[5].split("\t")[4]
            o.write(cluster[2][0:-1]+"**tsd1_"+tsd1+"_tsd2_"+tsd2+"_tir1_"+tir1+"_tir2_"+tir2+"\n")
            i=i+6
        i=i+1
    f.close()
    o.close()


getGff("%s_TIRvish.gff"%(genome_Name),"%s_TIRvish_pro.gff"%(genome_Name))

def GetFastaFromFile(argList):
    genomedb=argList[0]
    line=argList[1]
    outname=argList[2]
    gName=argList[3]
    o = open(outname,"a+") #shujun
    entry=line.split("\t")[0]
    p1 = int(line.split("\t")[3])
    p2 = int(line.split("\t")[4])
    infor=line.split("\t")[8].split("**")[-1]
    seq1 = subprocess.check_output("blastdbcmd -db '%s' -entry '%s' -range '%s'-'%s'" % (genomedb, entry, int(p1), int(p2)), shell=True)
    seq1 = seq1.decode("utf-8")
    out_seq1 = ''
    split = seq1.split('\n')
    for sp in split:
        if any([i.isdigit() for i in sp]):
            continue
        out_seq1 += sp
    o.write(">%s_%s_%s_%s_%s" % (gName,entry, p1, p2,infor) + str(out_seq1) + "\n") #shujun


if __name__ == '__main__':
    f="%s_TIRvish_pro.gff"%(genome_Name)
    outName=f[0:-4]+".fa"
    file=open(f,"r+")
    genomedb=genome_file+spliter+"db"
    genome_Name=genome_Name
    lines=file.readlines()
    argList=[[genomedb,i,outName,genome_Name] for i in lines]
    pool = multiprocessing.Pool(16)
    pool.map(GetFastaFromFile,argList)
    pool.close()
    pool.join()



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




def getSeq(file):
    remove=[]
    newf=open(file+spliter+"p","w")
    for rec in SeqIO.parse(file,"fasta"):
        s=str(rec.seq)
        #infor="_".join(rec.id.split("_")[4:])
        tirLen1=int(rec.id.split("_")[12])-int(rec.id.split("_")[11])+1
        tirLen2=int(rec.id.split("_")[15])-int(rec.id.split("_")[14])+1
        #tirLen = findDigitsSum(rec.id.split(":")[-2])
        tir1=str(rec.seq)[0:tirLen1]
        tir2=str(rec.seq)[0:tirLen2]
        if (TArepeats(s)!=True and checkNPer(s)!=True and TArepeats(tir1)!=True and checkN(tir1)!=True and TArepeats(tir2)!=True and checkN(tir2)!=True and  len(s)>=50 ):
            newf.write(">"+rec.description+"\n"+str(rec.seq)+"\n")
    newf.close()


file="%s_TIRvish_pro.fa"%(genome_Name)
getSeq(file)

mv = "mv *pro.fa ../temp/"
os.system(mv)
cp = "cp *-p ../temp/"
os.system(cp)





