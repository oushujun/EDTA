import pandas as pd
import os
import argparse


desired_width=500

pd.set_option('display.width', desired_width)

pd.set_option('display.max_columns',20)



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

targetDir=dir

spliter="-+-"


def split(s):
    return s.split("_")[1]+"_"+s.split("_")[2]+"_"+s.split("_")[3]


def FlankingBlast(file):
    f=pd.read_csv(file,header=None,sep="\t",comment="#")
   # f=pd.read_table(file,header=None,sep="\t",comment="#")
    f=f.loc[f[11]>=50]
    f["coor"]=f[0].apply(lambda x: split(x))
    coor=list(f["coor"])
    return list(set(coor))

def TArepeats(s):
    t=s.upper().count("T")
    a=s.upper().count("A")
    ta=t+a
    if (ta>len(s)*0.7):
        return True
    else:
        return False

def tir(s):
    return s.split()[0].split(":")[1]

def tsd(s):
    return s.split("_")[3].split(":")[1]

#os.chdir(targetDir)
os.chdir(targetDir+"/"+"%s"%(genome_Name)) #shujun

#for i in ["Module1","Module2","Module3"]:
#    cp= "cp %s/%s.gff3 ."%(genome_Name,genome_Name+spliter+i)
   # cp= "cp %s/%s/%s.gff3 ."%(i,genome_Name,genome_Name+spliter+i)
#    os.system(cp)
#    cp = "cp %s/%sLow ."%(genome_Name, genome_Name+spliter+i+spliter)
   # cp = "cp %s/%s/%sLow ."%(i, genome_Name, genome_Name+spliter+i+spliter)
#    os.system(cp)


############################################################################## Remove Entries with Homology in Flanking Sequences ##################################
print("############################################################################## Removing Entries with Homology in Flanking Sequences ##################################")

def removeIRFhomo(file,removelist,outputname):
    f=pd.read_table(file,header=None,sep="\t")
    print("%s entires in total"%(str(f.shape[0])))
    f["coor"]=f[0].apply(str)+"_"+f[3].apply(str)+"_"+f[4].apply(str)
    keep=f.loc[(~(f["coor"].isin(removelist)))]
    keep.to_csv(outputname,header=None,index=None,sep="\t")
    return keep

############################################################################## Remove Entries with Homology in Flanking Sequences ##################################

os.system("mkdir %s_combine"%(genome_Name))
for dataset in ["Module1","Module2","Module3"]:
    remove = FlankingBlast("%sLow"%(genome_Name+spliter+dataset+spliter))
    keep=removeIRFhomo("%s.gff3"%(genome_Name+spliter+dataset),remove,"%sClean.gff3"%(genome_Name+spliter+dataset+spliter))
    print("%s removed in %s" % (str(len(remove)), dataset))
    print("%s retained in %s"%(str(keep.shape[0]),dataset))
    cp ="cp %sClean.gff3 ./%s_combine"%(genome_Name+spliter+dataset+spliter,genome_Name)
    os.system(cp)
print("############################################################################## Finished Flanking Check ##################################")

######################################################################################################################################################################

def removeDupinSingle(file):
    f=pd.read_table(file,header=None,sep="\t")
    f=f.sort_values([0,3,4],ascending=[True,True,True])
    f=f.drop_duplicates([0,3,4],keep="last")
    return f

def RemoveTA(f):
    f["TIR"]=f[8].apply(lambda x:tir(x))
    f["TA"]=f["TIR"].apply(lambda x:TArepeats(x))
    sub=f.loc[f["TA"]==False]
    return sub[[0,1,2,3,4,5,6,7,8]]


def combineAll(f1,f2,f3,out):
    f12 = f1.append(f2, ignore_index=True)
    f123=f12.append(f3,ignore_index=True)
    f123.to_csv(out,header=None,index=None,sep="\t")
    return f123

############################################################################## Process and Combining three gff files ################################################
print("############################################################################## Process and Combining three gff files ################################################")

os.chdir("./%s_combine"%(genome_Name))
#
f_m1=removeDupinSingle("%sClean.gff3"%(genome_Name+spliter+"Module1"+spliter))
f_m2=removeDupinSingle("%sClean.gff3"%(genome_Name+spliter+"Module2"+spliter))
f_m3=removeDupinSingle("%sClean.gff3"%(genome_Name+spliter+"Module3"+spliter))
print("Removing TArepeats")
f_ta1=RemoveTA(f_m1)
f_ta2=RemoveTA(f_m2)
f_ta3=RemoveTA(f_m3)
print("Combine three modules")
comAll=combineAll(f_ta1, f_ta2,f_ta2, "%scombined.gff3"%(genome_Name+spliter))
print("#######################################################################Finished ##################################")
######################################################################################################################################################################

def ProcessGff(file,output):
    f=pd.read_table(file,header=None,sep="\t")
    f["pri"]=0
    mask=f[2]=="DTM"
    f.loc[mask, "pri"] = 1
    mask=f[2]=="DTC"
    f.loc[mask, "pri"] = 2
    mask = f[2] == "DTA"
    f.loc[mask, "pri"] = 3
    mask = f[2] == "DTT"
    f.loc[mask, "pri"] = 4
    mask = f[2] == "DTH"
    f.loc[mask, "pri"] = 5
    f["copy3"] = f[3]
    f["copy4"] = f[4]
    f.copy3 = f.copy3.shift(1).fillna(value=0).astype("int64")
    f.copy4 = f.copy4.shift(1).fillna(value=0).astype("int64")
    mask = ((f[3] == f["copy3"]) & (f[4] == f["copy4"]))
    f.loc[mask, 1] = "Both"
    f = f.sort_values([0, 3, 4, "pri",1], ascending=[True,True, True, True, True])
    f = f.drop_duplicates([0, 3, 4,"pri"], keep="first")
    f = f.sort_values([0, 3, 4,1, "pri"], ascending=[True,True,True, True, True])
    f=f.drop_duplicates([0,3,4],keep="first")
    f["length"]=f[4]-f[3]+1
    f=f[[0,1,2,3,4,5,6,7,8,"pri","length"]]
    f.to_csv(output,header=None,index=None,sep="\t")

############################################################################## Preparing for Removing Overlaps ################################################
print("############################################################################## Preparing for Removing Overlaps ##################################")

ProcessGff("%scombined.gff3"%(genome_Name+spliter), "%scombine_all_process.gff3"%(genome_Name+spliter))
f = pd.read_table("%scombine_all_process.gff3"%(genome_Name+spliter), header=None, sep="\t")
contig=list(set(list(f[0])))
for i in contig:
    chr=f.loc[f[0]==i]
    chr.to_csv("%s_combine_all_process_%s.gff3"%(genome_Name,i),header=None,index=None,sep="\t")
print("############################################################################## Finished: Removing Overlaps ##################################")
######################################################################################################################################################################

def splitInfor(s,i):
    l=s.split("_")
    return float(l[i])


def ProcessSelect(file):
    f = pd.read_table(file, header=None, sep="\t")
    f = f.sort_values([0,3, 4, 9,10], ascending=[True, True, True, True,True])

    f["copy3"]=f[3]
    f["copy4"]=f[4]
    f["copy10"]=f[10]
    f["TIRp"] = f[8].apply(lambda x: splitInfor(x, 2))
    f["TSDp"] = f[8].apply(lambda x: splitInfor(x, 5))
    f["copy9"] = f[9]
    f["copyTIRp"] = f["TIRp"]
    f["copyTSDp"] = f["TSDp"]
    f.copyTIRp = f.TIRp.shift(1).fillna(value=0)
    f.copyTSDp = f.TSDp.shift(1).fillna(value=0)
    f.copy3=f.copy3.shift(1).fillna(value=0).astype("int64")
    f.copy4 = f.copy4.shift(1).fillna(value=0).astype("int64")
    f.copy10 = f.copy10.shift(1).fillna(value=0).astype("int64")
    f.copy9 = f.copy9.shift(1).fillna(value=0).astype("int64")
    f.copyTIRp = f.TIRp.shift(1).fillna(value=0)
    f.copyTSDp = f.TSDp.shift(1).fillna(value=0)
    f["3_copy4"]=f[3]-f["copy4"]
    f["4_copy4"]=f[4]-f["copy4"]
    f["len-len"]=f[10]-f["copy10"]
    f["pri-pri"]=f[9]-f["copy9"]
    f["tir-tir"]=f["TIRp"]-f["copyTIRp"]
    f["tsd-tsd"]=f["TSDp"]-f["copyTSDp"]
    f=f.sort_values([0,3,4,1],ascending=[True,True,True,True])
    return f



def getRemoveList(f):
    removeList=[]
    overlap=(((f[3]-f["copy3"]<=30) & (f[3]-f["copy3"]>=0)) | ((f["copy4"]-f[4]<=30) & (f["copy4"]-f[4]>=0) ))
    r1=f.loc[(f["3_copy4"]<=0)&(f["4_copy4"]<=0) & (overlap==True)]
    if r1.shape[0]==0:
        pass
    removeList=removeList+list(r1.index.values)
    r2=f.loc[(f["3_copy4"]<=0)&(f["4_copy4"]>=0)]
    if r2.shape[0]==0:
        pass
    for index, row in r2.iterrows():
        if(row["pri-pri"]>0):
            removeList.append(index)
        elif(row["pri-pri"]<0):
            removeList.append(index-1)
        else:
            if(row["tir-tir"]<0):
                removeList.append(index)
            elif(row["tir-tir"]>0):
                removeList.append(index-1)
            else:
                if (row["tsd-tsd"] < 0):
                    removeList.append(index)
                elif (row["tsd-tsd"] > 0):
                    removeList.append(index-1)
                else:
                    if(row["len-len"] <= 0):
                        removeList.append(index)
                    else:
                        removeList.append(index - 1)
    return removeList


def CheckOverlap(file):
    re_open = pd.read_table(file, header=None, sep="\t")
    re_open["copy3"] = re_open[3]
    re_open["copy4"] = re_open[4]
    re_open["copy10"] = re_open[10]
    re_open.copy3 = re_open.copy3.shift(1).fillna(value=0).astype("int64")
    re_open.copy4 = re_open.copy4.shift(1).fillna(value=0).astype("int64")
    re_open.copy10 = re_open.copy10.shift(1).fillna(value=0).astype("int64")
    re_open["3_copy4"] = re_open[3] - re_open["copy4"]
    re_open["4_copy4"] = re_open[4] - re_open["copy4"]
    re_open["len-len"] = re_open[10] - re_open["copy10"]
    overlap = (((re_open[3]-re_open["copy3"]<=30) & (re_open[3]-re_open["copy3"]>=0)) | ((re_open["copy4"]-re_open[4]<=30) & (re_open["copy4"]-re_open[4]>=0) ))
    r1 = re_open.loc[(re_open["3_copy4"] <= 0) & (re_open["4_copy4"] <= 0) & (overlap==True)]
    r2 = re_open.loc[(re_open["3_copy4"] <= 0) & (re_open["4_copy4"] >= 0)]
    if (r1.shape[0]!=0 or r2.shape[0]!=0):
        return False
    else:
        return True

def deleteOverlap(file,output):
    f=ProcessSelect(file)
    l=getRemoveList(f)
    if len(l)!=0:
        newf = f.drop(f.index[l])
        newf = newf[[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]]
        newf.to_csv(output, header=None, index=None, sep="\t")
        deleteOverlap(output, output)
    else:
        newf = f[[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]]
        newf.to_csv(output, header=None, index=None, sep="\t")
        return newf

############################################################################## Removing Overlaps ################################################
print("############################################################################## Removing Overlaps ##################################")
f = pd.read_table("%scombined.gff3"%(genome_Name+spliter), header=None, sep="\t")
contig=list(set(list(f[0])))
for i in contig:
    newf = deleteOverlap("%s_combine_all_process_%s.gff3"%(genome_Name,i), "%s_combine_all_process_%s_Oremoved.gff3.txt" % (genome_Name,i))
    cat="cat *_combine_all_process_*_Oremoved.gff3.txt > %s_FinalAnn.gff3"%(genome_Name)
    os.system(cat)
print("############################################################################## Finished: Removing Overlaps ##################################")
######################################################################################################################################################################


def RemoveOerlap(file,outname):
    removeList=[]
    f=pd.read_table(file,header=None,sep="\t")
    f = f.sort_values([0, 3, 4, 9, 10], ascending=[True, True, True, True, True])

    f["copy3"] = f[3]
    f["copy4"] = f[4]
    f["copy10"] = f[10]
    f.copy3 = f.copy3.shift(1).fillna(value=0).astype("int64")
    f.copy4 = f.copy4.shift(1).fillna(value=0).astype("int64")
    f.copy10 = f.copy10.shift(1).fillna(value=0).astype("int64")
    f["3_copy4"] = f[3] - f["copy4"]
    f["4_copy4"] = f[4] - f["copy4"]
    f["len-len"] = f[10] - f["copy10"]
    f["copy_0"]=f[0]
    f["0_copy0"]=f[0]==f["copy_0"]

    f = f.sort_values([0, 3, 4, 1], ascending=[True, True, True, True])
    sub1=f.loc[(f["3_copy4"]<=0)&(f["len-len"]<=0)&(f["0_copy0"]==True)]
    removeList = removeList + list(sub1.index.values)
    sub2=f.loc[(f["3_copy4"]<=0)&(f["len-len"]>0)&(f["0_copy0"]==True)]
    removeList = removeList + [i-1 for i in list(sub2.index.values)]
    f=f.drop(f.index[removeList])
    f=f.drop(["copy3","copy4","copy10","3_copy4","4_copy4","len-len"],axis=1)
    f.to_csv(outname,header=None,index=None,sep="\t")
    return removeList

############################################################################## Deleting internal copies ################################################
#print("############################################################################## Deleting internal copies ##################################")

#RemoveOerlap("%s_FinalAnn.gff3"%(genome_Name), "%s_FinalAnn_Clint.gff3"%(genome_Name))

print("############################################################################## Finished: Deleting internal copies ##################################")
######################################################################################################################################################################
final=pd.read_table("%s_FinalAnn.gff3"%(genome_Name),header=None,sep="\t")
print(final[0:10])
final=final[list(range(0,9))]
final.to_csv("%s_FinalAnn.gff3"%(genome_Name),header=None,index=None,sep="\t")

#final=pd.read_table("%s_FinalAnn_Clint.gff3"%(genome_Name),header=0,sep="\t")
#final=final[list(range(0,9))]
#final.to_csv("%s_FinalAnn_Clint.gff3"%(genome_Name),header=None,index=None,sep="\t")

os.system("cp %s_FinalAnn.gff3 ../"%(genome_Name))
#os.system("cp %s_FinalAnn_Clint.gff3 ../"%(genome_Name))

print("############################################################################## Extracting Fasta ##################################")

