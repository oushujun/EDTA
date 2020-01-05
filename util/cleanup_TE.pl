#!/usr/bin/perl -w
use strict;
use FindBin;
use File::Basename;

# Clean up TE-related sequences from CDS
# Shujun Ou (shujun.ou.1@gmail.com) 10/28/2019

my $script_path = $FindBin::Bin;
my $cds = "";
my $rawlib = ""; #gene sequences in this file will be removed
my $minlen = 300; #minimal cds length to be retained
my $TEsorter = "$script_path/../bin/TEsorter/TEsorter.py";
my $output_by_list = "$script_path/output_by_list.pl";
my $cleanup = "$script_path/cleanup_tandem.pl";
my $name_code_decode = "$script_path/name_code_decode.pl";
my $threads = 4;
my $repeatmasker = "";

# read parameters
my $k=0;
foreach (@ARGV){
	$cds=$ARGV[$k+1] if /^-cds$/i and $ARGV[$k+1] !~ /^-/;
	$rawlib=$ARGV[$k+1] if /^-rawlib$/i and $ARGV[$k+1] !~ /^-/;
	$minlen=$ARGV[$k+1] if /^-minlen$/i and $ARGV[$k+1] !~ /^-/;
	$TEsorter=$ARGV[$k+1] if /^-tesorter$/i and $ARGV[$k+1] !~ /^-/;
	$threads=$ARGV[$k+1] if /^-threads$|^-t$/i and $ARGV[$k+1] !~ /^-/;
	$repeatmasker = $ARGV[$k+1] if /^-repeatmasker$/i and $ARGV[$k+1] !~ /^-/;
	$k++;
	}

# check files
die "The CDS file is empty or not exist!\n" unless -s $cds;
die "The raw library file is empty or not exist!\n" unless -s $rawlib;

# preprocess cds
my $cds_file = basename($cds);
`ln -s $cds ./` unless -e $cds_file;
$cds = $cds_file;
`perl -nle '\$_=(split)[0]; print \$_' $cds > $cds.mod`;
$cds = "$cds.mod";

# 1st attempt to find TEs in CDS with TEsorter
`python2 $TEsorter $cds -p $threads`;

# make an initial TE list
`cat $cds.rexdb.cls.tsv > $cds.TE.list`;
`grep -P "transposable|transposon|LINE" $cds >> $cds.TE.list`;

# get TE and non-TE seq from cds
`perl $output_by_list 1 $cds 1 $cds.TE.list -FA -ex > $cds.rmTE`;
`perl $output_by_list 1 $cds 1 $cds.TE.list -FA > $cds.TE`;

# convert seq names for RepeatMasker
`perl $name_code_decode 1 $cds.rmTE`;

# 2nd attempt to identify TEs in CDS based on repeatedness
`${repeatmasker}RepeatMasker -pa $threads -q -no_is -norna -nolow -div 40 -cutoff 225 -lib $cds.rmTE.code $rawlib 2>/dev/null`;
`awk '{print \$10}' $rawlib.out |sort|uniq -c|awk '{if (\$1>=10) print \$2}' | perl $output_by_list 1 $cds.rmTE.code 1 - -FA >> $cds.TE`; #CDS seqs appears >=10 times in masking the TE rawlib are considered TEs and removed from the CDS file
`awk '{print \$10}' $rawlib.out |sort|uniq -c|awk '{if (\$1>=10) print \$2}' | perl $output_by_list 1 $cds.rmTE.code 1 - -FA -ex > $cds.rmTE2`;

# 3rd attempt, mask remaining TE seqs in cds with identified TE seqs
if (-s "$cds.TE"){
	`${repeatmasker}RepeatMasker -pa $threads -q -no_is -norna -nolow -div 40 -lib $cds.TE -cutoff 225 $cds.rmTE2`;
	`perl $cleanup -Nscreen 1 -nc 300 -nc 0.3 -minlen $minlen -maxlen 300000 -cleanN 1 -cleanT 0 -trf 0 -f $cds.rmTE2.masked > $cds.noTE`;
	} else {
	print STDERR "\t\t\t\tWarning: No TE-related CDS found ($cds.TE empty). Will not use the self-cleaning step.\n\n";
	`cp $cds.rmTE.code $cds.noTE`;
	}

