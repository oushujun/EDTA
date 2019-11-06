#!/usr/bin/perl -w
use strict;
use FindBin;
use File::Basename;

# Clean up TE-related sequences from CDS
# Shujun Ou (shujun.ou.1@gmail.com) 10/28/2019

my $script_path = $FindBin::Bin;
my $cds = "";
my $minlen = 300; #minimal cds length to be retained
my $TEsorter = "$script_path/../bin/TEsorter/TEsorter.py";
my $output_by_list = "$script_path/output_by_list.pl";
my $cleanup = "$script_path/cleanup_tandem.pl";
my $name_code_decode = "$script_path/name_code_decode.pl";
my $threads = 4;

# read parameters
my $k=0;
foreach (@ARGV){
	$cds=$ARGV[$k+1] if /^-cds$/i and $ARGV[$k+1] !~ /^-/;
	$minlen=$ARGV[$k+1] if /^-minlen$/i and $ARGV[$k+1] !~ /^-/;
	$TEsorter=$ARGV[$k+1] if /^-tesorter$/i and $ARGV[$k+1] !~ /^-/;
	$threads=$ARGV[$k+1] if /^-threads$|^-t$/i and $ARGV[$k+1] !~ /^-/;
	$k++;
	}

# preprocess cds
my $cds_file = basename($cds);
`ln -s $cds ./` unless -e $cds_file;
$cds = $cds_file;
`perl -nle '\$_=(split)[0]; print \$_' $cds > $cds.mod`;
$cds = "$cds.mod";

# find TE-related seqs with TEsorter
`python2 $TEsorter $cds -p $threads`;

# make an initial TE list
`cat $cds.rexdb.cls.tsv > $cds.TE.list`;
`grep -P "transposable|transposon|LINE" $cds >> $cds.TE.list`;

# get TE and non-TE seq from cds
`perl $output_by_list 1 $cds 1 $cds.TE.list -FA -ex > $cds.rmTE`;
`perl $output_by_list 1 $cds 1 $cds.TE.list -FA > $cds.TE`;

# convert seq names for RepeatMasker
`perl $name_code_decode 1 $cds.rmTE`;

# mask remaining TE seqs in cds with identified TE seqs
if (-s "$cds.TE"){
	`RepeatMasker -pa $threads -q -no_is -norna -nolow -div 40 -lib $cds.TE -cutoff 225 $cds.rmTE.code`;
	`perl $cleanup -Nscreen 1 -nc 300 -nc 0.3 -minlen $minlen -maxlen 300000 -cleanN 1 -cleanT 0 -trf 0 -f $cds.rmTE.code.masked > $cds.noTE`;
	} else {
	print STDERR "\t\t\t\tWarning: No TE-related CDS found ($cds.TE empty). Will not use the self-cleaning step.\n\n";
	`cp $cds.rmTE.code $cds.noTE`;
	}


