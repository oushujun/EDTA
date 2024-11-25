#!/bin/bash -login
#SBATCH -N 1
#SBATCH -t 10-0:00:00
#SBATCH --ntasks-per-node 36
#SBATCH --mem=300GB

echo 'Start:'
date 

conda activate RM

path=/home/oushujun/las/TE/TEbenchmark/LTR/LtrDetector
genome=/home/oushujun/las/TE/TEbenchmark/Rice_MSU7.fasta
name=Rice_MSU7.fasta
rawlib=Rice_TIGR7Detector.bed.scn.list.fa

#goto local temp dir and make softlinks for files
cd $TMPDIR
ln -s $path/$rawlib $rawlib
ln -s $genome $name

#clean tandem
perl ~/las/git_bin/TElib_benchmark/util/cleanup_tandem.pl -misschar N -nc 50000 -nr 0.9 -minlen 80 -minscore 3000 -cleanN 1 -cleanT 1 -trf 1 -f $rawlib > $rawlib.noTandem
perl -nle 's/\|.*//; print $_' $rawlib.noTandem > $rawlib.noTandem.rename
rawlib=$rawlib.noTandem.rename

#clean nested
/usr/bin/time -v perl ~/las/git_bin/LTR_retriever/bin/cleanup_nestedIN.pl -in $rawlib -threads 30 > $rawlib.cln
/usr/bin/time -v perl ~/las/git_bin/LTR_retriever/bin/cleanup_nestedIN.pl -in $rawlib.cln -threads 30 > $rawlib.cln2
/usr/bin/time -v perl ~/las/git_bin/LTR_retriever/bin/cleanup_nestedIN.pl -in $rawlib.cln2 -threads 30 > $rawlib.cln3
/usr/bin/time -v perl ~/las/git_bin/LTR_retriever/bin/cleanup_nestedIN.pl -in $rawlib.cln3 -threads 30 > $rawlib.cln4
/usr/bin/time -v perl ~/las/git_bin/LTR_retriever/bin/cleanup_nestedIN.pl -in $rawlib.cln4 -threads 30 > $rawlib.cln5

lib=$rawlib.cln5

/usr/bin/time -v RepeatMasker -e ncbi -pa 36 -q -no_is -norna -nolow -div 40 -lib $lib -cutoff 225 $name
cp -r * $path
cd $path

/usr/bin/time -v perl ~/las/git_bin/TElib_benchmark/lib-test.pl -genome $name -std ../../Rice_MSU7.fasta.std6.9.2.out -tst $name.out -cat LTR -unknown 1

echo 'End:'
date

