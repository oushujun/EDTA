
# The Extensive *de novo* TE Annotator (EDTA)

## Introduction
This package is developed for automated whole-genome *de-novo* TE annotation and benchmarking the annotation performance of TE libraries.

For the initial search of TE candidates, [LTRharvest](http://genometools.org/), [LTR_FINDER_parallel](https://github.com/oushujun/LTR_FINDER_parallel), and [LTR_retriever](https://github.com/oushujun/LTR_retriever) are incorporated in this package to identify LTR retrotransposons; [TIR-Learner](https://github.com/weijiaweijia/TIR-Learner-Rice) and [MITE-Hunter](http://target.iplantcollaborative.org/mite_hunter.html) are incorporated in this package to identify TIR transposons (a subclass of DNA transposons); [HelitronScanner](https://sourceforge.net/projects/helitronscanner/) is incorporated in this package to identify *Helitron* transposons (a subclass of DNA transposons); and finally [RepeatModeler](https://github.com/rmhubley/RepeatModeler) is used to identify any TEs missed by these structure-based programs.

The EDTA package was designed to filter out false discoveries in raw TE candidates and generate a high-quality non-redundant TE library for whole-genome TE annotation.

For benchmarking of a testing TE library, I have provided the curated TE annotation (v6.9.5) for the rice genome (TIGR7/MSU7 version). You may use the `lib-test.pl` script to compare the annotation performance of your method/library to the methods we have tested (usage shown below).

## Installation
    conda create -n EDTA
    conda activate EDTA
    conda install -y -c conda-forge perl perl-text-soundex multiprocess regex
    conda install -y -c cyclus java-jdk
    conda install -y -c bioconda cd-hit repeatmodeler
    conda install -y -c bioconda/label/cf201901 repeatmasker
    conda install -y -c anaconda biopython pandas glob2
    conda install -y -c anaconda scikit-learn=0.19.0
    git clone https://github.com/oushujun/EDTA
    ./EDTA/EDTA.pl

## EDTA Usage
Activate the EDTA program:

    conda activate EDTA

Form head to toe (you got a genome and you want to get a high-quality TE library):
    
    perl EDTA.pl -genome your_genome.fasta -threads 36

Just the body (you got raw TE candidates from various programs and you want to filter them using EDTA):

    perl EDTA_process.pl [options]
      -genome	[File]	The genome FASTA
      -ltr	[File]	The raw LTR library FASTA
      -tir	[File]	The raw TIR library FASTA
      -mite	[File]	The raw MITE library FASTA
      -helitron	[File]	The raw Helitron library FASTA
      -repeatmasker [path]	The directory containing RepeatMasker (default: read from ENV)
      -blast [path]	The directory containing Blastn (default: read from ENV)
      -threads	[int]	Number of theads to run this script
      -help|-h	Display this help info

## Benchmarking
If you got a TE library and want to compare it's annotation performance to other methods, you can:

1.annotate the rice genome with your test library:

    RepeatMasker -pa 36 -q -no_is -norna -nolow -div 40 -lib custom.TE.lib.fasta -cutoff 225 rice_genome.fasta

2.Test the annotation performance of a particular TE category.

    perl lib-test.pl -genome genome.fasta -std genome.stdlib.RM.out -tst genome.testlib.RM.out -cat [options]
        -genome	[file]	FASTA format genome sequence
        -std	[file]	RepeatMasker .out file of the standard library
        -tst	[file]	RepeatMasker .out file of the test library
        -cat	[string]	Testing TE category. Use one of LTR|nonLTR|LINE|SINE|TIR|MITE|Helitron|Total|Classified
        -N	[0|1]	Include Ns in total length of the genome. Defaule: 0 (not include Ns).
        -unknown	[0|1]	Include unknown annotations to the testing category. This should be used when
                        the test library has no classification and you assume they all belong to the
                        target category specified by -cat. Default: 0 (not include unknowns)

eg.

    perl lib-test.pl -genome rice_genome.fasta -std ./EDTA/database/Rice_MSU7.fasta.std6.9.5.out -tst rice_genome.fasta.test.out -cat LTR


## Other resources
You may download the [rice genome here](http://rice.plantbiology.msu.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_7.0/all.dir/all.con).
