#!/usr/bin/env perl
use strict;
use warnings;
use FindBin;
use File::Basename;
# Shujun Ou (shujun.ou.1@gmail.com) 08/05/2020

my $usage = "\nMake a pangenome TE library from a list of individual libraries.
Files in -liblist ending with 'mod.EDTA.TElib.novel.fa' are treated as novel libraries.
Novel libraries will not be purged by -curatedlib as it is already novel.\n
perl make_panTElib.pl -liblist TElib.list [options]
	-liblist	[file]	A list of individual TE libraries (with accessible paths). One file per line.
	-curatedlib	[file]	A curated TE library. Optional. TEs from -liblist that are new to this lib will be kept.
	-miniden	[int]	Clustering parameter. Minimal identity in range (1, 100]. Default 80.
	-mincov		[float]	Clustering parameter. Minimal coverage in range [0, 1]. Default 0.95.
	-minlen		[int]	Clustering parameter. Minimal length in range [1, inf). Default 80.
	-threads	[int]	Number of threads. Default 4.
	-help			Display this info and quit.
";

# predefined parameters
my $liblist = '';
my $HQlib = '';
my $min_iden = "80"; #0-100
my $min_cov = "0.95"; #0-1
my $min_len = "80"; #>0
my $threads = 4;
my $debug = 0;

# dependencies
my $script_path = $FindBin::Bin;
my $cleanup_tandem = "$script_path/cleanup_tandem.pl";
my $cleanup_nested = "$script_path/cleanup_nested.pl";
my $rename_TE = "$script_path/rename_TE.pl";
my $repeatmasker = '';

my $k=0;
foreach (@ARGV){
	$liblist = $ARGV[$k+1] if /^-liblist$/i and $ARGV[$k+1] !~ /^-/;
	$HQlib = $ARGV[$k+1] if /^-curatedlib$/i and $ARGV[$k+1] !~ /^-/;
	$min_iden = $ARGV[$k+1] if /^-miniden$/i and $ARGV[$k+1] !~ /^-/;
	$min_cov = $ARGV[$k+1] if /^-mincov$/i and $ARGV[$k+1] !~ /^-/;
	$min_len = $ARGV[$k+1] if /^-minlen$/i and $ARGV[$k+1] !~ /^-/;
	$repeatmasker = $ARGV[$k+1] if /^-repeatmasker$/i and $ARGV[$k+1] !~ /^-/;
	$threads = $ARGV[$k+1] if /^-t$|^-threads$/i and $ARGV[$k+1] !~ /^-/;
	$debug = 1 if /^-debug$/i;
	die $usage if /^-h$|^-help$/i;
	$k++;
	}

die $usage unless -s $liblist;
die "Error: specified -mincov is not in range [0, 1]\n!" unless $min_cov =~ /^[0-9.]+$/ and $min_cov >= 0 and $min_cov <= 1;
die "Error: specified -miniden is not in range (1, 100]\n!" unless $min_iden =~ /^[0-9.]+$/ and $min_iden > 1 and $min_iden <= 100;

# define RepeatMasker -pa parameter
my $rm_threads = int($threads/4);

# print info
chomp (my $date = `date`);
print "\n$date\nStart to combine individual TE libraries from $liblist.\n";
print "Parameters to cluster sequences:\n\tMinimum identity: $min_iden\n\tMinimum coverage: $min_cov\n\tMinimum length: $min_len\n";
print "Curated TE library is specified: $HQlib.\n" if -s $HQlib;

open List, "<$liblist" or die $usage;
my $count = 0;
`rm $liblist.panTE.raw` if -s "$liblist.panTE.raw";
while (<List>){
	chomp (my $lib = $_);
	die "The library file $lib in $liblist is not found!\n" unless -s $lib;
	my $lib_file = basename($lib);
	`ln -s $lib $lib_file` unless -e $lib_file;
	print "Working on individual library: $lib_file\n";

	# get novel seqs based on -curatedlib
	if (-s $HQlib and $lib_file !~ /mod.EDTA.TElib.novel.fa/){
		print "\tIdentify novel sequences in this library.\t" if -s $HQlib;
		`${repeatmasker}RepeatMasker -e ncbi -pa $rm_threads -q -no_is -norna -nolow -div 40 -lib $HQlib $lib_file 2>/dev/null`;
		`perl $cleanup_tandem -misschar N -nc 50000 -nr $min_cov -minlen $min_len -minscore 3000 -trf 0 -cleanN 1 -cleanT 0 -f $lib_file.masked > $lib_file.novel.fa`;
		`rm $lib_file.masked $lib_file.cat.gz`;
		$lib_file = "$lib_file.novel.fa";
		print "Done.\n";
		}
	`perl $rename_TE $lib_file $count >> $liblist.panTE.raw`;
	chomp ($count += `grep -c \\> $lib_file`);
	$count += 100;
	}

# rename existing panTE lib files
if (-e "$liblist.panTE.fa"){
	my $old_lib = `ls -l $liblist.panTE.fa|perl -nle 'my (\$month, \$day, \$time) = (split)[6,7,8]; \$time =~ s/://; print "\${month}_\${day}_\$time"'`;
	chomp $old_lib;
	print "$liblist.panTE.fa exists in this folder, renamed to $liblist.panTE.fa_$old_lib\n";
	`mv $liblist.panTE.fa $liblist.panTE.fa_$old_lib`;
	}

# remove redundant
print "Cluster all library sequences...\n";
`perl $cleanup_nested -in $liblist.panTE.raw -cov $min_cov -minlen $min_len -miniden $min_iden -t $threads`;
`perl $rename_TE $liblist.panTE.raw.cln > $liblist.panTE.fa`;
`cat $HQlib >> $liblist.panTE.fa` if -s $HQlib;

# check
if (-s "$liblist.panTE.fa"){
	print "Done!\n";
	} else {
	print "Job failed!\n";
	}
chomp ($date = `date`);
print "$date\n";
