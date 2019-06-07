from Bio.Seq import Seq
from Bio import SeqIO
import os
import multiprocessing
from multiprocessing import Pool
import argparse

#
TSD = {}
TSD["DTA"] = [8]
TSD["DTC"] = [3,2]
TSD["DTH"] = [3]
TSD["DTM"] = [9,8,7]
TSD["DTT"] = [2]

#DTA:8
#DTC:2/3
#DTH:3(twa)
#DTM:7-10
#DTT:2(TA)

def compare(tir1, tir2):
    d = 0
    for i in range(0, len(tir1)):
        if (tir1[i] != tir2[i]):
            d += 1
    return d


def slidingWindow(seq1, seq2, tsdlength):
    set1 = []
    set2 = []
    for i in range(0, len(seq1) - tsdlength + 1):
        set1.append(seq1[i:i + tsdlength])
        set2.append(seq2[i:i + tsdlength])
    return set1, set2


def ConservedTIR(s1):
    if (s1.lower()[0:4] == "cact" or s1.lower()[0:4] == "gtga"):
        return True
    else:
        return False

def ConservedDTH(set1,tsd_dffset,l):
    for i in tsd_dffset:
        if (tsd_dffset[i]<l*0.2):
            s1=set1[int(i.split(":")[0])]
            if (s1.lower()=="tta" or s1.lower=="taa"):
                return True
    return False

def ConservedDTT(set1,tsd_dffset,l):
    for i in tsd_dffset:
        if (tsd_dffset[i]<l*0.2):
            s1=set1[int(i.split(":")[0])]
            if (s1[0:2].lower()=="ta"):
                return True
    return False

def GetDiff(set1, set2):
    tsd_diff = {}
    for i in range(0, len(set1)):
        for j in range(0, len(set2)):
            name = str(i) + ":" + str(j)
            diff = compare(set1[i], set2[j])
            tsd_diff[name] = diff
    return tsd_diff


def isTSD(tsd_dffset, l):
    for i in tsd_dffset:
        if tsd_dffset[i] < l * 0.2:
            return True
    return False


def CheckTIR(arglist):
    rec=arglist[0]
    family=arglist[1]
    dic = {}
    s = str(rec.seq)[200:-200]
    l_List = list(range(10, int(len(s) / 2)))
    for l in l_List:
        s1 = s[0:l]
        s2_ = s[-l:]
        s2 = Seq(s2_).reverse_complement()
        d = compare(s1, s2)
        if family=="DTC":
            if d < l * 0.2 and ConservedTIR(s1)==True:
                dic[str(rec.id)] = l
                break
        else:
            if d < l * 0.2:
                dic[str(rec.id)] = l
                break
    return dic


def CheckTSD(arglist):
    rec=arglist[0]
    family=arglist[1]
    dic = {}
    s = str(rec.seq)
    l = TSD[family]
    for i in l:
        s1 = s[200-i:200]
        last20 = s[-200:]
        s2 = last20[0:i]
        set1, set2 = slidingWindow(s1, s2, i)
        dff = GetDiff(set1, set2)
        if (family == "DTH" and ConservedDTH(set1, dff, i) == True):
            dic[rec.id] = i
            break
        elif (family == "DTT" and ConservedDTT(set1, dff, i) == True):
            dic[rec.id] = i
            break
        elif (family != "DTH" and family !="DTT"):
            TSDexist = isTSD(dff, i)
            if (TSDexist == True):
                dic[rec.id] = i
                break
    return dic

def getTSD(tsd_dffset, l_tsd, set1, set2):
    for i in tsd_dffset:
        if tsd_dffset[i] < l_tsd * 0.2:
            seq1 = set1[int(i.split(":")[0])]
            seq2 = set2[int(i.split(":")[1])]
            return seq1, seq2

def TIRpercent(seq1, seq2):
    d = compare(seq1, seq2)
    l = len(seq1)
    p = (l - d) / l
    p = p * 100
    p = round(p, 2)
    return p


def TSDpercent(seq1, seq2):
    d = compare(seq1, seq2)
    l = len(seq1)
    p = (l - d) / l
    p = p * 100
    p = round(p, 2)
    return p


def writeTofa(file, both, withTIR, withTSD):
    used = []
    output=file[0:-2]+"checkedM3"+".fa"
    w = open(output, "w")
    record = list(SeqIO.parse(file, "fasta"))
    for rec in record:
        if rec.id in both and rec.id not in used:
            s = str(rec.seq)[200:-200]
            l_tir = withTIR[rec.id]
            s1 = s[0:l_tir]
            s2_ = s[-l_tir:]
            s2 = Seq(s2_).reverse_complement()
            s2 = str(s2)
            p_tir = TIRpercent(s1, s2)
            s = str(rec.seq)
            l_tsd=withTSD[rec.id]
            s1tsd = s[200-l_tsd:200]
            last200 = s[-200:]
            s2tsd = last200[0:l_tsd]
            set1, set2 = slidingWindow(s1tsd, s2tsd, l_tsd)
            dff = GetDiff(set1, set2)
            seq1, seq2 = getTSD(dff, l_tsd, set1, set2)
            pTSD = TSDpercent(seq1, seq2)
            w.write(">" + str(rec.id) + "-+-"+"TIR:" + str(s1) + "_" + str(s2_) + "_" + str(p_tir) + "_" + "TSD:" + str(
                seq1) + "_" + str(seq2) + "_" + str(pTSD) + "\n" + str(rec.seq)[200:-200] + "\n")
            used.append(str(rec.id))

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

if __name__ == '__main__':
    os.chdir(targetDir)
    for i in ["DTA","DTC","DTH","DTM","DTT"]:
        cat="cat *%s > %s"%(i+".fa",genome_Name+spliter+i+spliter+"candidatesM3.fa")
        os.system(cat)
 #       rm ="rm *_%s.fa"%(i)
 #       os.system(rm)
    os.system("rm *.csv")
    os.system("rm *_nonHomo.fa")
    for i in ["DTA","DTC","DTH","DTM","DTT"]:
        file=genome_Name+spliter+i+spliter+"candidatesM3.fa"
        records=list(SeqIO.parse(file,"fasta"))
        records=[rec for rec in records if len(str(rec.seq))>=450]
        arglists=[[rec,file.split(spliter)[-2]] for rec in records]
        pool1 = multiprocessing.Pool(int(t))
        L_tir = pool1.map(CheckTIR, arglists)
        pool1.close()
        pool1.join()
        pool2 = multiprocessing.Pool(int(t))
        L_tsd = pool2.map(CheckTSD, arglists)
        pool1.close()
        pool1.join()
        withTIR={}
        for d in L_tir:
            if len(d)!=0:
                withTIR.update(d)
        withTSD={}
        for d in L_tsd:
            if len(d)!=0:
                withTSD.update(d)
        TIR_key=[k for k in withTIR]
        TSD_key=[k for k in withTSD]
        both = set(TIR_key).intersection(set(TSD_key))
        if (len(both)>0):
            writeTofa(file, both, withTIR, withTSD)
