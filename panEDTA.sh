#!/usr/bin/env bash

# This is the serial version of panEDTA
#	1. Each genome will be annotated by EDTA sequentially
#	2. Individual TE libraries will be combined by panEDTA
#	3. Each genome will be reannotated by the panEDTA library
# Shujun Ou (shujun.ou.1@gmail.com) 
# 06/21/2023 v0.2
# 10/10/2022 v0.1


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
   echo "\nPan-genome annnotation of TEs using EDTA"
   echo "Usage: bash $0 -g genome_list.txt -c cds.fasta -t 10"
   echo "   -g	A list of genome files with paths accessible from the working directory.
			Required: You can provide only a list of genomes in this file (one column, one genome each row).
			Option 1: You can also provide both genomes and CDS files in this file (two columns, one genome and 
				  one CDS each row). Missing of CDS files (eg, for some or all genomes) is allowed.
			Option 2: If there are existing EDTA annotations for individual genomes located in the directory of 
				  genome files, these annotations will be copied to the current work directory and skip 
				  recreating individual EDTA annotations."
   echo "   -c		Optional. Coding sequence files in fasta format.
   			The CDS file provided via this parameter will fill in the missing CDS files in the genome list.
			If no CDS files are provided in the genome list, then this CDS file will be used on all genomes."
   echo "   -l	Optional. A manually curated, non-redundant library following the RepeatMasker naming format."
   echo "   -f	Minimum number of full-length TE copies in individual genomes to be kept as candidate TEs for the pangenome.
   			Lower is more inclusive, and will ↑ library size, ↑ sensitivity, and ↑ inconsistency.
			Higher is more stringent, and will ↓ library size, ↓ sensitivity, and ↓ inconsistency.
			Default: 3."
   echo "   -a	Optional. Just generate the panEDTA library (0) or 
   			Perform whole-genome annotation using the generated panEDTA library (default, 1)."
   echo "   -o	Optional. Overwrite EDTA results. If you run panEDTA in the same folder containing EDTA results, the program
			will abort (default, 0) or overwrite (1)."
   echo "   -t	Number of CPUs to run panEDTA. Default: 10."
   echo ""
   exit 1 # Exit script after printing help
}

# Preset parameters
genome_list=''
cds=''
curatedlib=''
overwrite=0
fl_copy=3
threads=10
anno=1

# Read user inputs
while getopts "g:c::l::f::a::o::t::" opt
do
   case "$opt" in
      g ) genome_list="$OPTARG" ;;
      c ) cds="$OPTARG" ;;
      l ) curatedlib="$OPTARG" ;;
      f ) fl_copy="$OPTARG" ;;
      a ) anno="$OPTARG" ;;
      o ) overwrite="$OPTARG" ;;
      t ) threads="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Check parameters
if [ ! -s "$genome_list" ]; then
   echo "ERROR: The genomes $genome_list file is not found or is empty";
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
outfile=$(basename $genome_list 2>/dev/null)

### Begin panEDTA and print all parameters
printf "\n%s\n" "$(date)"
echo "Pan-genome Extensive de-novo TE Annotator (panEDTA)"
echo "   Output directory: $dir"
echo "   Genome files: $genome_list"
echo "   Coding sequences: $cds"
echo "   Curated library: $curatedlib"
echo "   Copy number cutoff: $fl_copy"
echo "   Overwrite EDTA results: $overwrite"
echo "   CPUs: $threads"
echo ""

## Step 1, initial EDTA annotation, consider to add --sensitive 1, consider to submit each EDTA job to different nodes.
# make softlink to global cds
if [ $cds != '' ]; then
	cds_file=$cds
	cds=`basename $cds_file 2>/dev/null`
	cds_file=`realpath $cds_file`
	ln -s $cds_file $cds 2>/dev/null
fi

# process one line each time
genomes="" # store a list of genomes
#cat "$genome_list" | while IFS= read -r i; do
while IFS= read -r i; do
	genome_file=`echo $i|awk '{print $1}'`

	# skip empty lines
	if [ "$genome_file" = "" ]; then
		break
	fi

	# check if genome file exist
	if [ ! -s "$genome_file" ]; then
		echo "ERROR: $genome_file specified by -g not exist!"
		exit 1
	fi

	# make softlink to genome
	genome_file=`realpath $genome_file`
	genome=`basename $genome_file`
	genomes="$genomes $genome"
	if [ ! -s $genome ]; then
		ln -s $genome_file $genome 2>/dev/null
	fi

	# make softlink to cds
	cds_ind_file=`echo $i|awk '{print $2}'`
	cds_ind=`basename $cds_ind_file 2>/dev/null` 
	if [ "$cds_ind" != '' ]; then
		if [ ! -s $cds_ind ]; then
			echo "ERROR: $cds_ind specified by -g not exist!"
			exit 1
		fi
		cds_ind_file=`realpath $cds_ind_file`
		ln -s $cds_ind_file $cds_ind 2>/dev/null
	fi

	# use the global $cds to replace a missing cds
	if [ "$cds_ind" = '' ] && [ "$cds" != '' ]; then
		cds_ind=$cds
	fi


	# check if current folder has EDTA results
	if [ `realpath "$genome_file.mod.EDTA.TEanno.sum"` = `realpath "$genome.mod.EDTA.TEanno.sum"` ] && [ $overwrite = 0 ]; then
		echo "ERROR: Existing EDTA result found for $genome and the Overwrite parameter (-o) is $overwrite!"
		exit 1
	fi

	# check if provided genome has EDTA annotation in the same folder
	if [ -s "$genome_file.mod.EDTA.TEanno.sum" ] && [ ! -s "$genome.mod.EDTA.TEanno.sum" ]; then # link annotations to the work directory
		echo "Existing EDTA annotation found in the directory of $genome, will use this as the panEDTA input"
		ln -s "$genome_file.mod" "$genome.mod" 2>/dev/null
		ln -s "$genome_file.mod.EDTA.TEanno.sum" "$genome.mod.EDTA.TEanno.sum" 2>/dev/null
		ln -s "$genome_file.mod.EDTA.TElib.fa" "$genome.mod.EDTA.TElib.fa" 2>/dev/null
		ln -s "$genome_file.mod.EDTA.TElib.novel.fa" "$genome.mod.EDTA.TElib.novel.fa" 2>/dev/null
		mkdir "$genome.mod.EDTA.raw" 2>/dev/null
		ln -s "$genome_file.mod.EDTA.raw/$genome.mod.RM2.fa" "$genome.mod.EDTA.raw/" 2>/dev/null
		ln -s "$genome_file.mod.EDTA.raw/$genome.mod.EDTA.intact.fa" "$genome.mod.EDTA.raw/" 2>/dev/null
		ln -s "$genome_file.mod.EDTA.raw/$genome.mod.EDTA.intact.gff3" "$genome.mod.EDTA.raw/" 2>/dev/null
		mkdir "$genome.mod.EDTA.combine" 2>/dev/null
		ln -s "$genome_file.mod.EDTA.combine/$genome.mod.EDTA.fa.stg1" "$genome.mod.EDTA.combine/" 2>/dev/null

		mkdir "$genome.mod.EDTA.anno" 2>/dev/null
		ln -s "$genome_file.mod.EDTA.anno/$genome.mod.EDTA.RM.out" "$genome.mod.EDTA.anno/" 2>/dev/null
	fi

#if [ ! 0 ]; then #test
	# de novo EDTA if no existing annotation is found
	if [ ! -s "$genome.mod.EDTA.TEanno.sum" ]; then
	# run a new EDTA if annotation of the genome is not existing
	echo "De novo annotate genome $genome with EDTA"
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
#fi #test
done < "$genome_list"

## Step 2, make pan-genome lib (quick step, use a sigle node is fine)
# get fl-TE with ≥ $fl_copy copies in each genome
printf "\n%s\n" "$(date)"
for genome in $genomes; do
	printf "\tIdenfity full-length TEs for genome %s\n"
	perl $path/util/find_flTE.pl $genome.mod.EDTA.anno/$genome.mod.EDTA.RM.out | \
		awk '{print $10}'| \
		sort| \
	 	uniq -c |\
		perl -snale 'my ($count, $id) = (split); next if $count < $fl_copy; print $_' -- -fl_copy=$fl_copy | \
		awk '{print $2"#"}' > $genome.mod.EDTA.TElib.fa.keep.list
done

# extract pan-TE library candidate sequences
printf "\n%s\n\tExtract pan-TE library candidate sequences\n" "$(date)"
for genome in $genomes; do
	if [ -s "$curatedlib" ]; then
		# a) if --curatedlib is provided
		for j in $(cat "$genome.mod.EDTA.TElib.fa.keep.list"); do
#		cat "$genome.mod.EDTA.TElib.fa.keep.list" | while read -r j; do
			grep "$j" "$genome.mod.EDTA.TElib.novel.fa"; 
		done | \
			perl "$path/util/output_by_list.pl" 1 "$genome.mod.EDTA.TElib.novel.fa" 1 - -FA > "$genome.mod.EDTA.TElib.fa.keep.ori"

#		while read -r j; do
#			grep "$j" "$genome.mod.EDTA.TElib.novel.fa";
#		done < "$genome.mod.EDTA.TElib.fa.keep.list" | \
#			perl "$path/util/output_by_list.pl" 1 "$genome.mod.EDTA.TElib.novel.fa" 1 - -FA > "$genome.mod.EDTA.TElib.fa.keep.ori"

	else
		# b) if --curatedlib is not provided
	#	for j in $(cat "$genome.mod.EDTA.TElib.fa.keep.list"); do 
		for j in `cat "$genome.mod.EDTA.TElib.fa.keep.list"`; do 
	#	cat "$genome.mod.EDTA.TElib.fa.keep.list" | while read -r j; do
			grep "$j" "$genome.mod.EDTA.TElib.fa"; 
		done | \
        		perl "$path/util/output_by_list.pl" 1 "$genome.mod.EDTA.TElib.fa" 1 - -FA > "$genome.mod.EDTA.TElib.fa.keep.ori"

#		while read -r j; do
#			grep "$j" "$genome.mod.EDTA.TElib.fa";
#		done < "$genome.mod.EDTA.TElib.fa.keep.list" | \
#			perl "$path/util/output_by_list.pl" 1 "$genome.mod.EDTA.TElib.fa" 1 - -FA > "$genome.mod.EDTA.TElib.fa.keep.ori"

	fi
done

# aggregate TE libs
i=0
for genome in $genomes; do
	i=$(($i+5000)); 
	perl $path/util/rename_TE.pl $genome.mod.EDTA.TElib.fa.keep.ori $i; 
done | perl $path/util/rename_TE.pl - > $outfile.panEDTA.TElib.fa.raw

# remove redundant
printf "\n%s\n\tGenerate the panEDTA library\n" "$(date)"
perl $path/util/cleanup_nested.pl -in $outfile.panEDTA.TElib.fa.raw -cov 0.95 -minlen 80 -miniden 80 -t $threads
perl -nle 's/>(TE_[0-9]+)/>pan$1/; print $_' $outfile.panEDTA.TElib.fa.raw.cln > $outfile.panEDTA.TElib.fa

# Extra step if --curatedlib is provided:
if [ -s "$curatedlib" ]; then
	cat $curatedlib >> $outfile.panEDTA.TElib.fa
fi
printf "\n%s\n" "$(date)"
printf "\tpanEDTA library of %s is generated: %s.panEDTA.TElib.fa\n" "$genome_list" "$outfile"
#printf "\tPan-genome library: %s.panEDTA.TElib.fa\n" "$outfile"

## Step 3, re-annotate all genomes with the panEDTA library, consider to submit each RepeatMasker and EDTA job to different nodes.
if [ "$anno" = '1' ]; then
	for genome in $genomes; do
		printf "\n%s\nReannotate genome %s with the panEDTA library - homology\n" "$(date)" $genome
		#echo "Reannotate genome $genome with the panEDTA library - homology"
		if [ ! -s "$genome.mod.panEDTA.out" ]; then
			ln -s $genome.mod $genome.mod.panEDTA
			RepeatMasker -e ncbi -pa $rm_threads -q -div 40 -lib $genome_list.panEDTA.TElib.fa -cutoff 225 -gff $genome.mod.panEDTA >/dev/null
		fi
		perl -i -nle 's/\s+DNA\s+/\tDNA\/unknown\t/; print $_' $genome.mod.panEDTA.out
	done

while IFS= read -r i; do
	genome=`basename $(echo $i|awk '{print $1}') 2>/dev/null`
	cds_ind=`basename $(echo $i|awk '{print $2}') 2>/dev/null`

	# skip empty lines
        if [ $genome = '' ]; then
                break
        fi

        # use the global $cds to replace a missing cds
        if [ "$cds_ind" = '' ] && [ "$cds" != '' ]; then
                cds_ind=$cds
        fi

	echo "Reannotate genome $genome with the panEDTA library - structural"
	perl $path/EDTA.pl --genome $genome -t $threads --step final --anno 1 --curatedlib $genome_list.panEDTA.TElib.fa --cds $cds_ind --rmout $genome.mod.panEDTA.out
done < $genome_list

printf "%s\n\tpanEDTA annotation of %s is finished!\n\n" "$(date)" "$genome_list"
fi
