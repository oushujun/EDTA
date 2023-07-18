#!/usr/bin/env bash

# This is the serial version of panEDTA
#	1. Each genome will be annotated by EDTA sequentially
#	2. Individual TE libraries will be combined by panEDTA
#	3. Each genome will be reannotated by the panEDTA library
# Shujun Ou (shujun.ou.1@gmail.com) 
# 10/10/2022


## Test files explained
# Col.test.fa	Genome 1
# Ler.test.fa	Genome 2
# genome.cds.fa	CDS sequence. Could be absent, or from all input genomes, or from one of the genomes, or from a closely related species.
#		If users provided only one CDS file, then use it for all genomes. Each genome could use a different CDS.
#		This file is optional, but having it will improve the annotation by purging gene sequences from TEs.
# athrep.ref.formatted	Existing TE library from this species. Optional. Providing a curated library will significantly improve
#		the annotation and also harvest existing knowledge.


# Help info
helpFunction()
{
   echo ""
   echo "Usage: $0 -g genome_list.txt -c cds.fasta -t 10"
   echo -e "\t-g	A list of genome files with paths accessible from the working directory.
			Required: You can provide only a list of genomes in this file (one column, one genome each row).
			Optional: You can also provide both genomes and CDS files in this file (two columns, one genome and one CDS each row).
				  Missing of CDS files (eg, for some or all genomes) is allowed."
   echo -e "\t-c		Optional. Coding sequence files in fasta format.
   			The CDS file provided via this parameter will fill in the missing CDS files in the genome list.
			If no CDS files are provided in the genome list, then this CDS file will be used on all genomes."
   echo -e "\t-l	Optional. A manually curated, non-redundant library following the RepeatMasker naming format."
   echo -e "\t-f	Minimum number of full-length TE copies in individual genomes to be kept as candidate TEs for the pangenome.
   			Lower is more inclusive, and will ↑ library size, ↑ sensitivity, and ↑ inconsistency.
			Higher is more stringent, and will ↓ library size, ↓ sensitivity, and ↓ inconsistency.
			Default: 3."
   echo -e "\t-a	Optional. Just generate the panEDTA library (0) or 
   			Perform whole-genome annotation using the generated panEDTA library (default, 1)."
   echo -e "\t-t	Number of CPUs to run panEDTA. Default: 10."
   echo ""
   exit 1 # Exit script after printing help
}

# Preset parameters
genome_list=''
cds=''
curatedlib=''
fl_copy=3
threads=10
anno=1

# Read user inputs
while getopts "g:c::l::f::a::t::" opt
do
   case "$opt" in
      g ) genome_list="$OPTARG" ;;
      c ) cds="$OPTARG" ;;
      l ) curatedlib="$OPTARG" ;;
      f ) fl_copy="$OPTARG" ;;
      a ) anno="$OPTARG" ;;
      t ) threads="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Check parameters
if [ ! -s "$genome_list" ]; then
   echo "ERROR: The genomes $genomes_list file is not found or is empty";
   helpFunction
fi

if [ "$cds" != '' ] && [ ! -s "$cds" ]; then
   echo "ERROR: The cds $cds file is not found or is empty"
   helpFunction
fi

if [ "$curatedlib" != '' ] && [ ! -s "$curatedlib" ]; then
   echo "ERROR: The curated library $curatedlib file is not found or is empty"
   helpFunction
fi

# Set paths
path=$(dirname "$0") #program path
dir=$(pwd) #current work directory
rm_threads=$(($threads/4))

### Begin panEDTA and print all parameters
echo `date`
echo "Pan-genome Extensive de-novo TE Annotator (panEDTA)"
echo -e "\tGenome files: $genome_list"
echo -e "\tCoding sequences: $cds"
echo -e "\tCurated library: $curatedlib"
echo -e "\tCopy cutoff: $fl_copy"
echo -e "\tCPUs: $threads"

## Step 1, initial EDTA annotation, consider to add --sensitive 1, consider to submit each EDTA job to different nodes.
# make softlink to global cds
if [ $cds != '' ]; then
	cds_file=$cds
	cds=`basename $cds_file 2>/dev/null`
	ln -s $cds_file $cds 2>/dev/null
fi

# process one line each time
IFS="
"
for i in `cat $genome_list`; do
	genome_file=`echo $i|awk '{print $1}'`

	# skip empty lines
	if [ $genome_file == '' ]; then
		break
	fi

	# make softlink to genome
	genome=`basename $genome_file`
	if [ ! -s $genome ]; then
		ln -s $genome_file $genome 2>/dev/null
	fi

	# make softlink to cds
	cds_ind_file=`echo $i|awk '{print $2}'`
	cds_ind=`basename $cds_ind_file 2>/dev/null` 
	if [ ! -s $cds_ind ] && [ $cds_ind != '' ]; then
		ln -s $cds_ind_file $cds_ind 2>/dev/null
	fi

	# use the global $cds to replace a missing cds
	if [ "$cds_ind" == '' ] && [ "$cds" != '' ]; then
		cds_ind=$cds
	fi

	echo "Annotate genome $genome with EDTA"
	if [ ! -s "$genome.mod.EDTA.TEanno.sum" ]; then # run a new EDTA if annotation of the genome is not existing
		perl $path/util/count_base.pl $genome -s > $genome.stats
		if [ "$curatedlib" != '' ]; then
			perl $path/EDTA.pl --genome $genome -t $threads --anno 1 --cds $cds_ind --curatedlib $curatedlib || {
				echo >&2 ERROR: Initial EDTA failed for $genome
			exit 1
		}
		else
			perl $path/EDTA.pl --genome $genome -t $threads --anno 1 --cds $cds_ind || {
				echo >&2 ERROR: Initial EDTA failed for $genome
			exit 1
		}
		fi
	fi
done

## Step 2, make pan-genome lib (quick step, use a sigle node is fine)
# get fl-TE with ≥ $fl_copy copies in each genome
for i in `cat $genome_list`; do   
	genome=`basename $(echo $i|awk '{print $1}') 2>/dev/null`

	# skip empty lines
        if [ $genome == '' ]; then
                break
        fi

	echo "Idenfity full-length TEs for genome $genome"
	perl $path/util/find_flTE.pl $genome.mod.EDTA.anno/$genome.mod.EDTA.RM.out | \
		awk '{print $10}'| \
		sort| \
	 	uniq -c |\
		perl -snale 'my ($count, $id) = (split); next if $count < $fl_copy; print $_' -- -fl_copy=$fl_copy | \
		awk '{print $2"#"}' > $genome.mod.EDTA.TElib.fa.keep.list
done

# extract pan-TE library candidate sequences
for i in `cat $genome_list`; do
	genome=`basename $(echo $i|awk '{print $1}') 2>/dev/null`

	# skip empty lines
        if [ $genome == '' ]; then
                break
        fi

	if [ -s "$curatedlib" ]; then

		# a) if --curatedlib is provided
		for j in `cat $genome.mod.EDTA.TElib.fa.keep.list`; do
			grep $j $genome.mod.EDTA.TElib.novel.fa; 
		done | \
			perl $path/util/output_by_list.pl 1 $genome.mod.EDTA.TElib.novel.fa 1 - -FA > $genome.mod.EDTA.TElib.fa.keep.ori
	else
		# b) if --curatedlib is not provided
		for j in `cat $genome.mod.EDTA.TElib.fa.keep.list`; do 
			grep $j $genome.mod.EDTA.TElib.fa; 
		done | \
        		perl $path/util/output_by_list.pl 1 $genome.mod.EDTA.TElib.fa 1 - -FA > $genome.mod.EDTA.TElib.fa.keep.ori
	fi
done

# aggregate TE libs
i=0
for j in `cat $genome_list`; do
	genome=`basename $(echo $j|awk '{print $1}') 2>/dev/null`

	# skip empty lines
        if [ $genome == '' ]; then
                break
        fi

	i=$(($i+5000)); 
	perl $path/util/rename_TE.pl $genome.mod.EDTA.TElib.fa.keep.ori $i; 
done | perl $path/util/rename_TE.pl - > $genome_list.panEDTA.TElib.fa.raw

# remove redundant
echo "Generate the panEDTA library"
perl $path/util/cleanup_nested.pl -in $genome_list.panEDTA.TElib.fa.raw -cov 0.95 -minlen 80 -miniden 80 -t $threads
cp $genome_list.panEDTA.TElib.fa.raw.cln $genome_list.panEDTA.TElib.fa

# Extra step if --curatedlib is provided:
if [ -s "$curatedlib" ]; then
	cat $curatedlib >> $genome_list.panEDTA.TElib.fa
fi
echo `date`
echo -e "\tpanEDTA library of $genome_list is generated!"

## Step 3, re-annotate all genomes with the panEDTA library, consider to submit each RepeatMasker and EDTA job to different nodes.
if [ "$anno" == '1' ]; then
for i in `cat $genome_list`; do
	genome=`basename $(echo $i|awk '{print $1}') 2>/dev/null`

	# skip empty lines
        if [ $genome == '' ]; then
                break
        fi

	echo "Reannotate genome $genome with the panEDTA library - homology"
	RepeatMasker -pa $rm_threads -q -div 40 -lib $genome_list.panEDTA.TElib.fa -cutoff 225 -gff $genome.mod >/dev/null
	perl -i -nle 's/\s+DNA\s+/\tDNA\/unknown\t/; print $_' $genome.mod.out
done

for i in `cat $genome_list`; do
	genome=`basename $(echo $i|awk '{print $1}') 2>/dev/null`
	cds_ind=`basename $(echo $i|awk '{print $2}') 2>/dev/null`

	# skip empty lines
        if [ $genome == '' ]; then
                break
        fi

        # use the global $cds to replace a missing cds
        if [ "$cds_ind" == '' ] && [ "$cds" != '' ]; then
                cds_ind=$cds
        fi

	echo "Reannotate genome $genome with the panEDTA library - structural"
	perl $path/EDTA.pl --genome $genome -t $threads --step final --anno 1 --curatedlib $genome_list.panEDTA.TElib.fa --cds $cds_ind --rmout $genome.mod.out
done

echo `date`
echo -e "\tpanEDTA annotation of $genome_list is finished!"
fi
