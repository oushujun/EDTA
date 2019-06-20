from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os
import multiprocessing
from multiprocessing import Pool
import argparse
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


def GetFastaFromFile(argsList):
    line=argsList[0]
    genomeName=argsList[1]
    outname=genomeName+spliter+"F200.fa"
    genomeFile=argsList[2]
    p1 = int(line.split("\t")[3])
    p2 = int(line.split("\t")[4])
    entry=line.split("\t")[0]
    family=str(line.split("\t")[2])
    o=open(outname, 'a+')
    if (p1 < p2):
        p_start = p1-200
        p_end = p2+200
    else:
        p_start = p2-200
        p_end = p1+200
    seq1 = subprocess.check_output("blastdbcmd -db '%s' -entry '%s' -range '%s'-'%s'" % (genomeFile+spliter+"db", entry, int(p_start), int(p_end)), shell=True)

    seq1 = seq1.decode("utf-8")
    out_seq1 = ''
    split = seq1.split('\n')
    for sp in split:
        if any([i.isdigit() for i in sp]):
            continue
        out_seq1 += sp
    o.write(">"+genomeName+spliter+entry+spliter+str(p1)+spliter+str(p2)+spliter+str(p_start)+spliter+str(p_end)+spliter+family + "\n" + str(out_seq1) + "\n")


if __name__ == '__main__':
    f=open(genome_Name+spliter+"Module1.gff3","r+")
    lines=f.readlines()
    l=[[i,genome_Name,genome_file] for i in lines]
    pool = multiprocessing.Pool(int(t))
    pool.map(GetFastaFromFile,l)
    pool.close()
    pool.join()


def getLTR(rec):
    s = str(rec.seq)
    ID = str(rec.id)
    seq1 = s[0:200]
    seq2 = s[-200:]
    seq1 = SeqRecord(Seq(seq1), id="seq1" + ID)
    seq2 = SeqRecord(Seq(seq2), id="seq2" + ID)
    SeqIO.write(seq1, "seq1_" + ID + ".fasta", "fasta")
    SeqIO.write(seq2, "seq2_" + ID + ".fasta", "fasta")
    blast = "blastn -query %s -subject %s -outfmt '7 qacc sacc length pident gaps mismatch qstart qend sstart send evalue qcovhsp' -out %s" % ("seq1_" + ID + ".fasta", "seq2_" + ID + ".fasta", ID + "_blast.tsv") #shujun
    os.system(blast)
#    rm = "rm %s" % ("seq2_" + ID + ".fasta")
 #   os.system(rm)
  #  rm = "rm %s" % ("seq1_" + ID + ".fasta")
   # os.system(rm)

if __name__ == '__main__':
    records = list(SeqIO.parse(genome_Name+spliter+"F200.fa", "fasta"))
    pool = multiprocessing.Pool(int(t))
    pool.map(getLTR, records)
    pool.close()
    pool.join()

def Remove(file):
    f=open(file,"r+")
    lines=f.readlines()
    l=[]
    for line in lines:
        if (line[0]!="#"):
            l.append(line)
    if len(l)==0:
       rm="rm %s 2>/dev/null"%(file)
       os.system(rm)

if __name__ == '__main__':

    files=os.listdir(".")
    filelist=[i for i in files if i[-4:]==".tsv"]
    pool = multiprocessing.Pool(int(t))
    pool.map(Remove,filelist)
    pool.close()
    pool.join()

cat = "cat *.tsv > %sModule1%sLow"%(genome_Name+spliter,spliter)
os.system(cat)
rm = "rm *.tsv seq1_* seq2_* 2>/dev/null" #shujun
os.system(rm)
