import pandas as pd
import os
import argparse



parser = argparse.ArgumentParser()#pylint: disable=invalid-name
parser.add_argument("-name", "--genomeName", help="Genome Name", required=True)
parser.add_argument("-p", "--path", help="Source code path", required=True)
parser.add_argument("-d", "--currentD", help="Path of current directory", required=True)
args = parser.parse_args()#pylint: disable=invalid-name


genome_Name = args.genomeName
path=args.path
dir=args.currentD

targetDir=dir+"/"+genome_Name+"/"
spliter="-+-"

def getFile(genomeName):
    for i in ["DTA", "DTC", "DTH", "DTM", "DTT"]:
        oldname=targetDir+genomeName+spliter+"blast"+spliter+"Rice_%s_RefLib" % (i)
        newname=oldname+spliter+"tem"
        f = open(oldname , "r+")
        o = open(newname , "w")
        lines = f.readlines()
        for line in lines:
            if (line[0] != "#"):
                o.write(line)
        f.close()
        o.close()
        mv = "mv %s %s"%(newname,oldname)
        os.system(mv)

getFile(genome_Name)

def ProcessHomology(genomeName):
    for i in ["DTA", "DTC", "DTH", "DTM", "DTT"]:
        blast = targetDir+ genomeName+spliter+"blast"+spliter+"Rice_%s_RefLib" % (i)
        f = pd.read_csv(blast, header=None, sep="\t") #shujun
      #  f = pd.read_table(blast, header=None, sep="\t")
        contigName=list(set(list(f[1])))
        for contig in contigName:
            con_f=f.loc[f[1]==contig]
            con_f = con_f.loc[(f[11] == 100) & (f[3] >= 80)]
            con_f = con_f.sort_values([1, 8, 9, 11, 3], ascending=[True, True, True, True, True])
            con_f = con_f.drop_duplicates([1, 8, 9], keep="last")
            con_f.to_csv(targetDir+contig+spliter+i+spliter+"select.csv", header=None, index=None, sep="\t")


ProcessHomology(genome_Name)

names=targetDir+genome_Name+spliter+"blast"+spliter+"Rice_%s_RefLib" % ("*")
#rm = "rm %s"%(names) #shujun
#os.system(rm) #shujun




