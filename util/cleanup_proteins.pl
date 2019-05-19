#!/usr/bin/perl -w
use strict;
use File::Basename;

## Shujun Ou (shujun.ou.1@gmail.com, 05/13/2019, Department of EEOB, ISU)

my $usage = "
	Clean up DNA TE and LINE related protein-coding sequences and Plant protein-coding sequences.
	Dependency: blastx is required accessible through \$PATH
	Usage: perl cleanup_proteins.pl sequence.fa
	Customization: You may change parameters by editing this file.
	\n";
	

my $target = $ARGV[0];
die $usage unless -s $target;

my $blastplus = '';
my $threads = 36;
my $script_path = "~/las/git_bin/LTR_retriever";
my $LINE="$script_path/database/Tpases020812LINE";
my $DNA="$script_path/database/Tpases020812DNA";
my $PlantP="$script_path/database/alluniRefprexp082813";
my $procovTE=0.7; #for DNA TE and LINE alignment, propotional, cumulated alignment coverage more than this will be treated as protein contained
my $procovPL=0.7; #for plant protein alignment, propotional, cumulated alignment coverage more than this will be treated as protein contained
my $prolensig=30; #aa, hits with alignment length less than this number will not be counted.

#prepare blastx databases
my $rand=int(rand(1000000));
my ($LINE_base, $LINE_path)=fileparse($LINE);
my ($DNA_base, $DNA_path)=fileparse($DNA);
my ($PlantP_base, $PlantP_path)=fileparse($PlantP);
`cp $LINE ./$LINE_base.$rand`;
`cp $DNA ./$DNA_base.$rand`;
`cp $PlantP ./$PlantP_base.$rand`;
$LINE="$LINE_base.$rand";
$DNA="$DNA_base.$rand";
$PlantP="$PlantP_base.$rand";

`${blastplus}makeblastdb -in $LINE -dbtype prot`;
`${blastplus}makeblastdb -in $DNA -dbtype prot`;
`${blastplus}makeblastdb -in $PlantP -dbtype prot`;


##DNA TE transposase, LINE transposase, and plant protein masking
`${blastplus}blastx -word_size 3 -outfmt 6 -max_target_seqs 10 -num_threads $threads -query $target -db $LINE -out $target.line.out`;
`${blastplus}blastx -word_size 3 -outfmt 6 -max_target_seqs 10 -num_threads $threads -query $target -db $DNA -out $target.dna.out`;
`cat $target.line.out $target.dna.out > $target.otherTE.out`;
`perl $script_path/bin/purger.pl -blast $target.otherTE.out -seq $target -cov $procovTE -purge 0 -len $prolensig`;

`${blastplus}blastx -word_size 3 -outfmt 6 -max_target_seqs 10 -num_threads $threads -query $target.clean -db $PlantP -out $target.plantP.out`;
`perl $script_path/bin/purger.pl -blast $target.plantP.out -seq $target.clean -cov $procovPL -purge 1 -len $prolensig`;

#clean up
`rm $LINE* $DNA* $PlantP*`;

