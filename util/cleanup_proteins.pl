#!/usr/bin/perl -w
use strict;
use File::Basename;
use FindBin;

## Shujun Ou (shujun.ou.1@gmail.com, 05/13/2019, Department of EEOB, ISU)

my $usage = "
	Clean up DNA TE and LINE related protein-coding sequences and/or Plant protein-coding sequences.
	Dependency: blastx is required accessible through \$PATH
	Usage: perl cleanup_proteins.pl -seq sequence.fa [options]
		-rmdnate	[0|1]	Clean up DNA TE transposase sequences (default: 0, not remove)
		-rmline		[0|1]	Clean up LINE transposase sequences (default: 0, not remove)
		-rmprot		[0|1]	Clean up protein-coding sequences (default: 1, remove plant coding proteins).
						For other systems (such as archaea, bacteria, fungi, invertebrates, 
						vertebrates, viruses), you must provide your own -Protlib.
		-protlib	[file]	Protein-coding aa sequences to be removed from -seq. (default lib: alluniRefprexp082813 (plant))
						You may use uniprot_sprot database available from here: 
						ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/
		-blast		[Path]	Path to the BLAST+ folder (Default: obtain from ENV)
		-purger		[File]	The purge.pl scrit to remove BLAST-matching regions (default: found in the same folder of this script).
		-threads	[int]	The number of threads to run this script. (default, 4)
	For more customizations, you may change internal parameters by editing this file.
	\n";

# pre-defined parameters
my $name = '';
my $target = '';
my $rmDNATE = 0;
my $rmLINE = 0;
my $rmProt = 1;

my $script_path = $FindBin::Bin;
my $LINE="$script_path/../database/Tpases020812LINE";
my $DNA="$script_path/../database/Tpases020812DNA";
my $Protlib="$script_path/../database/alluniRefprexp082813";
my $blastplus = '';
my $purger = "$script_path/purger.pl";

my $procovTE=0.7; #for DNA TE and LINE alignment, propotional, cumulated alignment coverage more than this will be treated as protein contained
my $procovPL=0.7; #for plant protein alignment, propotional, cumulated alignment coverage more than this will be treated as protein contained
my $prolensig=30; #aa, hits with alignment length less than this number will not be counted.
my $threads = 4;

# read custom parameters
my $k=0;
foreach (@ARGV){
	$target = $ARGV[$k+1] if /^-seq$/i and $ARGV[$k+1] !~ /^-/;
	$rmDNATE = $ARGV[$k+1] if /^-rmdnate$/i and $ARGV[$k+1] !~ /^-/;
	$rmLINE = $ARGV[$k+1] if /^-rmline$/i and $ARGV[$k+1] !~ /^-/;
	$rmProt = $ARGV[$k+1] if /^-rmprot$/i and $ARGV[$k+1] !~ /^-/;
	$Protlib = $ARGV[$k+1] if /^-protlib$/i and $ARGV[$k+1] !~ /^-/;
	$blastplus = $ARGV[$k+1] if /^-blast$/i and $ARGV[$k+1] !~ /^-/;
	$threads = $ARGV[$k+1] if /^-threads$/i and $ARGV[$k+1] !~ /^-/;
	$k++;
	}

die "Warning: Nothing will be cleaned based on your parameters.\n" if $rmDNATE eq 0 and $rmLINE eq 0 and $rmProt eq 0;
die "Sequence file not found!\n$usage" unless -s $target;
$name = $target;

#prepare blastx databases
my $rand=int(rand(1000000));

# clean up DNA TE
if ($rmDNATE eq 1){
	my ($DNA_base, $DNA_path)=fileparse($DNA);
	`cp $DNA ./$DNA_base.$rand`;
	$DNA="$DNA_base.$rand";
	`${blastplus}makeblastdb -in $DNA -dbtype prot`;
	`${blastplus}blastx -word_size 3 -outfmt 6 -max_target_seqs 10 -num_threads $threads -query $target -db $DNA -out $target.dnate.out`;
	`perl $purger -blast $target.dnate.out -seq $target -cov $procovTE -purge 0 -len $prolensig`;
	`cp $target.clean $target.dnate_clean`;
	$target = "$target.dnate_clean";
	`rm $DNA*`;
	}

# clean up LINE
if ($rmLINE eq 1){
	my ($LINE_base, $LINE_path)=fileparse($LINE);
	`cp $LINE ./$LINE_base.$rand`;
	$LINE="$LINE_base.$rand";
	`${blastplus}makeblastdb -in $LINE -dbtype prot`;
	`${blastplus}blastx -word_size 3 -outfmt 6 -max_target_seqs 10 -num_threads $threads -query $target -db $LINE -out $target.line.out`;
	`perl $purger -blast $target.line.out -seq $target -cov $procovTE -purge 0 -len $prolensig`;
	`cp $target.clean $target.line_clean`;
	$target = "$target.line_clean";
	`rm $LINE*`;
	}

# clean up Proteins
if ($rmProt eq 1){
	my ($Prot_base, $Prot_path)=fileparse($Protlib);
	`cp $Protlib ./$Prot_base.$rand`;
	$Protlib="$Prot_base.$rand";
	`${blastplus}makeblastdb -in $Protlib -dbtype prot`;
	`${blastplus}blastx -word_size 3 -outfmt 6 -max_target_seqs 10 -num_threads $threads -query $target -db $Protlib -out $target.prot.out`;
	`perl $purger -blast $target.prot.out -seq $target -cov $procovPL -purge 1 -len $prolensig`;
	`cp $target.clean $target.prot_clean`;
	$target = "$target.prot_clean";
	`rm $Protlib*`;
	}

# rename the final clean sequence
`cp $target $name.clean`;

