from Bio import SeqIO
import argparse


parser = argparse.ArgumentParser()#pylint: disable=invalid-name
parser.add_argument("-g", "--genomeFile", help="Genome file in fasta format", required=True)
parser.add_argument("-name", "--genomeName", help="Genome Name", required=True)

args = parser.parse_args()#pylint: disable=invalid-name

genome_file = args.genomeFile
genome_Name = args.genomeName


def getGenomeInformation(genome):
    records=SeqIO.parse(genome,"fasta")
    names=[rec.id for rec in records]
    for name in names:
        if "-+-" in names:
            print("sequence name has special character")
        c=names.count(name)
        if c>1:
            print(name +" presents"+str(c)+" times in the genome file")
    
    print(str(len(names)) +" sequences in genome file %s"%(genome_Name))
    
getGenomeInformation(genome_file)
