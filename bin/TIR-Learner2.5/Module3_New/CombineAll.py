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

os.chdir(targetDir)

def split(s):
#    return s.split("_")[1]+"_"+s.split("_")[2]+"_"+s.split("_")[3]
    return s.split(spliter)[1]+"_"+s.split(spliter)[2]+"_"+s.split(spliter)[3]

def TArepeats(s):
    t=s.upper().count("T")
    a=s.upper().count("A")
    ta=t+a
    if (ta>=len(s)*0.7):
        return True
    else:
        return False

def tir(s):
    return s.split()[0].split(":")[1]

def tsd(s):
    return s.split("_")[3].split(":")[1]

os.chdir(targetDir)

for i in ["Module3"]:
    cp= "cp %s/%s.gff3 ."%(genome_Name,genome_Name+spliter+i)
    os.system(cp)

#print("######################################### Removing Entries with Homology in Flanking Sequences ######################################")

def removeDupinSingle(file):
    f=pd.read_csv(file,header=None,sep="\t") #shujun
    f=f.sort_values([0,3,4,1,2],ascending=[True,True,True,True,True])
    f=f.drop_duplicates([0,3,4],keep="first")
    f=f.drop_duplicates([0,3],keep="first")
    f=f.drop_duplicates([0,4],keep="first")
    return f

def RemoveTA(f):
    f["TIR"]=f[8].apply(lambda x:tir(x))
    f["TA"]=f["TIR"].apply(lambda x:TArepeats(x))
    sub=f.loc[f["TA"]==False]
    return sub[[0,1,2,3,4,5,6,7,8]]



print("############################################################ Process and Combining three gff files ##################################")

#
f_m3=removeDupinSingle("%s.gff3"%(genome_Name+spliter+"Module3"))
print("Removing TArepeats")
f_ta3=RemoveTA(f_m3)
f_ta3.to_csv("%s_FinalAnn.gff3"%(genome_Name),header=None,index=None,sep="\t")


def getEntry(coor,entry):
    if (entry=="contig"):
       return coor.split(spliter)[0]
    if (entry==1):
       return int(coor.split(spliter)[1])
    else:
       return int(coor.split(spliter)[2])


def setClustters(CoorList):
    this=0
    New_list=[]
    while this < len(CoorList)-1:
          next=this+1
          this_coor=CoorList[this]
          cluster=[this_coor]
          while (next<len(CoorList)):
              next_coor = CoorList[next]
              next_p1 = getEntry(next_coor, 1)
              this_p2 = getEntry(this_coor, 2)
              if (next_p1<=this_p2):
                  cluster.append(next_coor)
                  next+=1
              elif(next!=len(CoorList)-1):
                  break
              else:
                  New_list.append([CoorList[next]])
                  break
          this=next
          New_list.append(cluster)
    return New_list

def removeCoor(clusterList):
    RemainList=[]
    for cluster in clusterList:
        if len(cluster)==1:
            RemainList.append(cluster[0])
        else:
            l=[abs(getEntry(c,2)-getEntry(c,1)) for c in cluster]
            for c2 in cluster:
                if (abs(getEntry(c2,2)-getEntry(c2,1))==max(l)):
                    RemainList.append(c2)
                    break
    return RemainList

     
def RemoveOverlap(finalFile):
    f=pd.read_table(finalFile,header=None,sep="\t")
    f=f.sort_values([0,3,4,1,2])
    contigs=list(set(list(f[0])))
    f["coor"]=f[0].apply(str)+spliter+f[3].apply(str)+spliter+f[4].apply(str)
    newDF=pd.DataFrame()
    for cn in contigs:
        sub=f.loc[f[0]==cn]
        CoorList=list(sub["coor"])
        clusterList=setClustters(CoorList)
        remain=removeCoor(clusterList)
        filter_df=sub.loc[sub["coor"].isin(remain)]
        newDF=newDF.append(filter_df)
    newDF.to_csv("%s_FinalAnn_filter.gff3"%(genome_Name),index=None,header=None,sep="\t") 

finalFile="%s_FinalAnn.gff3"%(genome_Name)
RemoveOverlap(finalFile)


print("############################################################ Finished ###############################################################")


