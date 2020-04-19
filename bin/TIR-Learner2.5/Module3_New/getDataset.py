import multiprocessing
from multiprocessing import Pool
import argparse
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' #mute all tensorflow info, warnings, and error msgs. #shujun
#os.environ["KMP_AFFINITY"] = "compact" #mute all OpenMP warnings. #shujun
os.environ["KMP_WARNINGS"] = '0' #mute all OpenMP warnings. #shujun
import warnings
warnings.filterwarnings('ignore',category=FutureWarning) #mute tensorflow warnings #shujun
from Bio import SeqIO
import numpy as np
from tensorflow.python.keras.utils import to_categorical
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
from tensorflow.python.keras.models import load_model
import subprocess
import tensorflow as tf
tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)
tf.logging.set_verbosity(tf.logging.ERROR)

parser = argparse.ArgumentParser()#pylint: disable=invalid-name
parser.add_argument("-name", "--genomeName", help="Genome Name", required=True)
parser.add_argument("-p", "--path", help="Source code path", required=True)
parser.add_argument("-t", "--processer", help="Number of processer", required=True)
parser.add_argument("-d", "--currentD", help="Path of current directory", required=True)
parser.add_argument("-g", "--genomeFile",help="Path to the genome file", required=True)
args = parser.parse_args()#pylint: disable=invalid-name

genome_Name = args.genomeName
path=args.path
t=args.processer
dir=args.currentD
genomeFile=args.genomeFile

targetDir=dir+"/"+genome_Name+"/"
os.chdir(targetDir)
spliter="-+-"


def getSeqFragment(arglist):
    file=arglist[0]
    output=arglist[1]
    featureSize=arglist[2]
    records=list(SeqIO.parse(file,"fasta"))
    if (len(records)>0): 
        f=open(output,"w")   
        for rec in records:
            sequence=str(rec.seq)
            if (len(sequence)>=featureSize*2):
                s=sequence[0:featureSize]+sequence[-featureSize:]
            else:
                s1=sequence[0:int(len(sequence)/2)]
                s2=sequence[int(len(sequence)/2):]
                n1="N"*(featureSize-len(s1))
                n2="N"*(featureSize-len(s2))
                s1=s1+n1
                s2=n2+s2
                s=s1+s2 
            if all(i in ["A","T","C","G","N"] for i in list(s)):
                f.write(">"+rec.description+"\n"+s+"\n")
        f.close()



if __name__ == '__main__':
    files=os.listdir(".")
    files=[i for i in files if i.split(spliter)[-1]=="p"]
    output=[file+spliter+"toPre.fa" for file in files]
    featureSize=200
    lists=[[files[i],output[i],200] for i in range(0,len(files))]
    pool = multiprocessing.Pool(int(t))
    d = pool.map(getSeqFragment,lists)
    pool.close()
    pool.join()


integer_encoder = LabelEncoder()
one_hot_encoder = OneHotEncoder()
input_features = []


def getData(file):
    feature_integer_encoder = LabelEncoder()
    input_features = []
    records=SeqIO.parse(file,"fasta")
    l_seq=[str(rec.seq) for rec in records]
    records=SeqIO.parse(file,"fasta")
    l_target=[rec1.id.split("_")[-1] for rec1 in records]
    voc=["A","C","G","T","N"]
    feature_integer_encoder.fit(voc)
    sequences = list(filter(None, l_seq))
    for sequence in sequences:
        integer_encoded = feature_integer_encoder.transform(list(sequence))
        integer_encoded = np.array(integer_encoded).reshape(-1, 1)
        s = to_categorical(integer_encoded,num_classes=len(voc))
        input_features.append(s)
    inputfeatures=np.array(input_features)
    np.save(file+spliter+"features.npy",inputfeatures)


if __name__ == '__main__':
    files=os.listdir(".")
    files=[i for i in files if i.split(spliter)[-1]=="toPre.fa"]
    pool = multiprocessing.Pool(int(t))
    d = pool.map(getData,files)
    pool.close()
    pool.join()


def Predict(file):
    model = load_model(path+'/CNN0912.h5')
    npData=file+spliter+"features.npy"
    
    prefeature=np.load(npData)
    if (prefeature.shape[0]>=1):
        predicted_labels = model.predict(np.stack(prefeature))
        l_class=["DTA","DTC","DTH","DTM","DTT","NonTIR"]
        y_classes = predicted_labels.argmax(axis=-1)
        target_integer_encoder = LabelEncoder()
        target_integer_encoder.fit(l_class)
        target_integer_encoded = target_integer_encoder.transform(l_class)
        d={}
        for i in range(0,len(target_integer_encoded)):
            d[target_integer_encoded[i]]=l_class[i]
        preFile=open(file+spliter+"predi.fa","w")
        records=list(SeqIO.parse(file,"fasta"))
        y_name=[d[i] for i in y_classes]
        for i in range (0,len(records)):
            preFile.write(">"+records[i].id+"_"+y_name[i]+"\n"+str(records[i].seq)+"\n")
        preFile.close()
    else:
        print(file+" has no candidate to be classified")

    

if __name__ == '__main__':
   files=os.listdir(".")
   files=[i for i in files if i.split(spliter)[-1]=="toPre.fa"]
   pool = multiprocessing.Pool(int(t))
   d = pool.map(Predict,files)
   pool.close()
   pool.join()


def GetListFromFile(file):
    listName=file+spliter+"200.list" #shujun
    records=list(SeqIO.parse(file,"fasta"))
    lines=[rec.id for rec in records]
    o = open(listName,"a+")
    for infor in lines: #shujun
        entry=infor.split(":")[0]
        p1 = int(infor.split(":")[1])
        p2 = int(infor.split(":")[2])
        fam=infor.split("_")[-1]
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
        o.write(genome_Name+spliter+entry+spliter+str(p1)+spliter+str(p2)+spliter+str(p_start)+spliter+str(p_end)+"_"+fam+"_200" + "\t" + entry + ":" + str(p_start) + ".." + str(p_end) + "\n")
    o.close() #shujun


def GetFastaFromList(argList): #shujun
    genomeFile=argList[0]
    listName=argList[1]
    outName=argList[2]
    get_seq = "perl %s/call_seq_by_list2.pl %s -C %s -header 1 -out %s" % (path, listName, genomeFile, outName)
    subprocess.run(['/bin/bash', '-c', get_seq])
    clean_head = "perl -i -nle 's/^>.*\|/>/; print $_' %s" % (outName)
    subprocess.run(['/bin/bash', '-c', clean_head])


if __name__ == '__main__':
    files=os.listdir(".")
    prefiles=[i for i in files if i.split(spliter)[-1]=="predi.fa"]
    argList=[[genomeFile,prefiles[i]+spliter+"200.list",prefiles[i]+spliter+"200"] for i in range(0,len(prefiles))] #shujun
    pool = multiprocessing.Pool(int(t))
    pool.map(GetListFromFile,prefiles) #shujun
    pool.map(GetFastaFromList,argList) #shujun
    pool.close()
    pool.join()

