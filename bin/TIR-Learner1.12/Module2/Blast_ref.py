import os
import argparse
import multiprocessing
from multiprocessing import Pool
import pandas as pd

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
spliter="-+-"
targetDir=dir+"/"+genome_Name+"/"
os.chdir(targetDir)


def BlastRef(file):
    ref_list=["Rice_DTA_RefLib","Rice_DTC_RefLib","Rice_DTH_RefLib","Rice_DTM_RefLib","Rice_DTT_RefLib"]
    for refLib in ref_list:
        blast = "blastn -query %s -subject %s -outfmt '7 qseqid sseqid length pident gaps mismatch qstart qend sstart send evalue qcovhsp' -out %s" % (file , path+"/RefLib/"+refLib, targetDir+file+spliter+"blast"+spliter+refLib)
       # blast = "blastn -query %s -db %s -num_threads %s -outfmt '7 qacc sacc length pident gaps mismatch qstart qend sstart send evalue qcovhsp' -out %s" % (file , targetDir+refLib+spliter+"db", int(t), targetDir+file+spliter+"blast"+spliter+refLib) #shujun
        os.system(blast)



def ProcessBlast(file):
    f=open(file,"r+")
    lines=f.readlines()
    lines=[i for i in lines if i[0]!="#"]
    if len(lines)==0:
        rm="rm %s 2>/dev/null"%(file)
        os.system(rm)
    else:
        newlines=[i.split("\t") for i in lines]
        newf=pd.DataFrame(newlines)
        newf.to_csv(file,header=None,index=None,sep="\t")
        read=pd.read_csv(file,header=None,sep="\t") #shujun
      #  read=pd.read_table(file,header=None,sep="\t")
        read=read.loc[(read[11]>=80) & (read[3] >=80) ]
        read=read.sort_values([0,11,3],ascending=[True,True,True])
        read = read.drop_duplicates([0], keep="last")
        read.to_csv("%s80" % (file+spliter), header=None, index=None, sep="\t")
        rm="rm %s 2>/dev/null"%(file)
        os.system(rm)


if __name__ == '__main__':
    os.chdir(targetDir)
    files=os.listdir(".")
    files=[i for i in files if i.split(spliter)[-1]=="p"]
    pool = multiprocessing.Pool(int(t))
    pool.map(BlastRef,files)
    pool.close()
    pool.join()

    files=os.listdir(".")
#    files=[i for i in files if i[-7:]=="_RefLib"]
    files=[i for i in files if i.split(spliter)[0]==genome_Name and i[-7:]=="_RefLib"] #shujun
    pool = multiprocessing.Pool(int(t))
    pool.map(ProcessBlast,files)
    pool.close()
    pool.join()
    rm = "rm ./*GRFmite.fa 2>/dev/null"
    os.system(rm)

