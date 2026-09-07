#!/usr/bin/env perl
use warnings;
use strict;
use FindBin;
use File::Basename;

# Clean up TE-related sequences from CDS
# Shujun Ou (shujun.ou.1@gmail.com) 10/28/2019

my $script_path = $FindBin::Bin;
my $cds = "";
my $rawlib = ""; #sequences in this file will be use to remove TEs in $cds
my $minlen = 300; #minimal cds length to be retained
my $output_by_list = "$script_path/output_by_list.pl";
my $cleanup = "$script_path/cleanup_tandem.pl";
my $name_code_decode = "$script_path/name_code_decode.pl";
my $threads = 4;
my $TEsorter = "";
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

# check dependencies
foreach my $script ($output_by_list, $cleanup, $name_code_decode){
	die "ERROR: The helper script $script is not found!\n" unless -s $script;
	}
$TEsorter=`command -v TEsorter 2>/dev/null` if $TEsorter eq '';
$TEsorter=~s/TEsorter\n$//;
die "ERROR: TEsorter is not found in the path $TEsorter!\n" unless -X "${TEsorter}TEsorter";
$repeatmasker=`command -v RepeatMasker 2>/dev/null` if $repeatmasker eq '';
$repeatmasker=~s/RepeatMasker\n$//;
die "ERROR: RepeatMasker is not found in the path $repeatmasker!\n" unless -X "${repeatmasker}RepeatMasker";

# preprocess cds
my $cds_file = basename($cds);
&run_cmd("ln -s $cds ./") unless -e $cds_file;
$cds = $cds_file;
&run_cmd("perl $name_code_decode 1 $cds 2>&1");
$cds = "$cds.code";

# 1st attempt to find TEs in CDS with TEsorter
&run_cmd("${TEsorter}TEsorter $cds -p $threads 2>&1");

# make an initial TE list
&run_cmd("cat $cds.rexdb.cls.tsv 2>&1 > $cds.TE.list");
&run_cmd("grep -P \"transposable|transposon|LINE\" $cds 2>&1 >> $cds.TE.list", 1); #exit 1 (no match) is not an error

# get TE and non-TE seq from cds
&run_cmd("perl $output_by_list 1 $cds 1 $cds.TE.list -FA -ex 2>&1 > $cds.rmTE");
&run_cmd("perl $output_by_list 1 $cds 1 $cds.TE.list -FA 2>&1 > $cds.TE");

# 2nd attempt to identify TEs in CDS based on repeatedness
&run_cmd("${repeatmasker}RepeatMasker -e ncbi -pa $threads -q -no_is -nolow -div 40 -cutoff 225 -lib $cds.rmTE $rawlib 2>&1");
&run_cmd("awk '{print \$10}' $rawlib.out |sort|uniq -c|awk '{if (\$1>=10) print \$2}' | perl $output_by_list 1 $cds.rmTE 1 - -FA 2>&1 >> $cds.TE"); #CDS seqs appears >=10 times in masking the TE rawlib are considered TEs and removed from the CDS file
&run_cmd("awk '{print \$10}' $rawlib.out |sort|uniq -c|awk '{if (\$1>=10) print \$2}' | perl $output_by_list 1 $cds.rmTE 1 - -FA -ex 2>&1 > $cds.rmTE2");

# 3rd attempt, mask remaining TE seqs in cds with potential TE seqs identified in cds ($cds.TE)
if (-s "$cds.TE"){
	&run_cmd("${repeatmasker}RepeatMasker -e ncbi -pa $threads -q -no_is -nolow -div 40 -lib $cds.TE -cutoff 225 $cds.rmTE2 2>&1");
	if (-s "$cds.rmTE2.masked"){
		&run_cmd("perl $cleanup -Nscreen 1 -nc 300 -nc 0.3 -minlen $minlen -maxlen 300000 -cleanN 1 -cleanT 0 -trf 0 -f $cds.rmTE2.masked 2>&1 > $cds.noTE");
		} else {
		&run_cmd("cp $cds.rmTE2 $cds.noTE 2>&1");
		}
	} else {
	print STDERR "\t\t\t\tWarning: No TE-related CDS found ($cds.TE empty). Will not use the self-cleaning step.\n\n";
	&run_cmd("cp $cds.rmTE $cds.noTE 2>&1");
	}


# run an external command; die with the captured output unless it exits 0 (or with one of the
# benign exit codes in $ok, e.g. 1 for grep's "no match"). Note "cmd 2>&1 > file": stderr is
# captured while stdout still goes to the file.
sub run_cmd {
	my ($cmd, $ok) = @_;
	$ok = '' unless defined $ok;
	my $out = `$cmd`;
	my $rc = $? >> 8;
	if ($? != 0 and " $ok " !~ / $rc /){
		die "ERROR: command failed (exit $rc): $cmd\n$out\n";
		}
	return $out;
	}

