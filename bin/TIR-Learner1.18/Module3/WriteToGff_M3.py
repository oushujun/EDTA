from Bio import SeqIO
import os
import argparse


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
spliter="-+-"

def writeTogff(fastafile,Family,output,source):
    record=list(SeqIO.parse(fastafile,"fasta"))
    out=open(output,"a+")
    for rec in record:
        ID=str(rec.id).split(spliter)
        seq_ID=ID[1]
        p1=int(ID[2])
        p2=int(ID[3])
        type = Family
        if (p1 < p2 ):
            start = p1
            end = p2
            strand = "."
        else:
            start = p2
            end = p1
            strand = "."
        length = end - start + 1
        if length < 50:
            pass
        attribut = ID[-1] + spliter + str(length)
        out.write(seq_ID + "\t" + source + "\t" + type + "\t" + str(start) + "\t" + str(
            end) + "\t" + "." + "\t" + strand + "\t" + "." + "\t" + attribut + "\n")

    out.close()



os.chdir(targetDir)
files=os.listdir(".")
files=[i for i in files if i.split(spliter)[-1]=="candidatesM3.checkedM3.fa"]
for file in files:
    Family=file.split(spliter)[-2]
    output=file[0:-3]+".gff3"
    writeTogff(file, Family, output,"Module3")
cat = "cat *%s*%s*.gff3 > %sModule3.gff3"%(spliter,spliter,genome_Name+spliter)
os.system(cat)
os.system("rm *%s*%s*.gff3"%(spliter,spliter))


