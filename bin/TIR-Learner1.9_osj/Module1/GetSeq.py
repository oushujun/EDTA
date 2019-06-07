import subprocess
import multiprocessing
from multiprocessing import Pool
import argparse
import os


parser = argparse.ArgumentParser()#pylint: disable=invalid-name
parser.add_argument("-g", "--genomeFile", help="Genome file in fasta format", required=True)
parser.add_argument("-name", "--genomeName", help="Genome Name", required=True)
parser.add_argument("-p", "--path", help="Source code path", required=True)
parser.add_argument("-d", "--currentD", help="Path of current directory", required=True)
parser.add_argument("-t", "--processer", help="Number of processer", required=True)
args = parser.parse_args()#pylint: disable=invalid-name

genome_file = args.genomeFile
genome_Name = args.genomeName
path=args.path
t=args.processer
dir=args.currentD

targetDir=dir+"/"+genome_Name+"/"
spliter="-+-"


#mkDB="makeblastdb -in %s -out %s -parse_seqids -dbtype nucl"%(genome_file,genome_file+spliter+"db") #shujun
#os.system(mkDB) #shujun

def GetFastaFromFile(argList):
    genomedb=argList[0]
    selectfile=argList[1]
    f=open(selectfile,"r+")
    lines=f.readlines()
    if len(lines)>0:
        outname=targetDir+selectfile[:-4]+".fa"
        o = open(outname,"w")
        for line in lines:
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
            seq1 = subprocess.check_output("blastdbcmd -db '%s' -entry '%s' -range '%s'-'%s'" % (genomedb, entry, int(p_start), int(p_end)), shell=True)
            seq1 = seq1.decode("utf-8")
            out_seq1 = ''
            split = seq1.split('\n')
            for sp in split:
                if any([i.isdigit() for i in sp]):
                    continue
                out_seq1 += sp
            o.write(">"+genome_Name+spliter+entry+spliter+str(p1)+spliter+str(p2)+spliter+str(p_start)+spliter+str(p_end) + "\n" + str(out_seq1) + "\n")
        f.close()
        o.close()

# #
# #
if __name__ == '__main__':
    os.chdir(targetDir)
    fileList=os.listdir(".")
    fileList=[i for i in fileList if i[-3:]=="csv"]
    genomedb=genome_file+spliter+"db"
    argList=[[genomedb,selectfile] for selectfile in fileList]
    pool = multiprocessing.Pool(int(t))
    pool.map(GetFastaFromFile,argList)
    pool.close()
    pool.join()

