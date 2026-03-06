#!/usr/bin/env perl
use warnings;
use strict;
use FindBin;

#Shujun Ou (shujun.ou.1@gmail.com)
#v0.1 02/14/2020
#v0.2 03/05/2026 Add mincov and divergence titration

my $usage = "\nEvaluate annotation consistency for given annotations.
	perl evaluation.pl -genome genome.fa -anno RepeatMasker.out -maxcount [int] -mincov [float] -threads [int]
		-anno	The whole-genome TE annotation file in RepeatMasker .out format
		-maxcount	The maximum number of stat lines to obtain. Default: 100000. 0 = no limit
		-mincov	Minimum reciprocal coverage for the redun category. Default: 0.95
\n";

my $script_path = $FindBin::Bin;
my $genome = '';
my $RMout = '';
my $blast = '';
my $threads = 4;
my $maxcount = 100000;
my $mincov = 0.95;
my $call_seq = "$script_path/call_seq_by_list.pl";
my $cleanup_nested = "$script_path/cleanup_nested.pl";
my $count_nested = "$script_path/count_nested.pl";

my $k=0;
foreach (@ARGV){
        $genome=$ARGV[$k+1] if /^-genome$/i and $ARGV[$k+1] !~ /^-/;
        $RMout=$ARGV[$k+1] if /^-anno$/i and $ARGV[$k+1] !~ /^-/;
	$maxcount=$ARGV[$k+1] if /^-maxcount$/i and $ARGV[$k+1] !~ /^-/;
	$mincov=$ARGV[$k+1] if /^-mincov$/i and $ARGV[$k+1] !~ /^-/;
        $blast=$ARGV[$k+1] if /^-blast$/i and defined $ARGV[$k+1] and $ARGV[$k+1] !~ /^-/;
        $threads=$ARGV[$k+1] if /^-threads$|^-t$/i and $ARGV[$k+1] !~ /^-/;
        $k++;
        }

#check files
die $usage unless -s $genome;
die $usage unless -s $RMout;

#get blast path
$blast=`command -v blastn 2>/dev/null` if $blast eq '';
$blast=~s/blastn\n//i;
$blast="$blast/" if $blast ne '' and $blast !~ /\/$/;
die "blastn is not exist in the BLAST+ path $blast!\n" unless -X "${blast}blastn";

my $date=`date`;
chomp ($date);
print "$date\tEvaluation starts...\n";

# extract whole-genome TE and perform all-v-all blast, then summarize the results
`awk '{if (\$5~/[0-9]+/ && \$1>300 && \$7-\$6>80) print \$11"\t"\$5":"\$6".."\$7}' $RMout | perl $call_seq - -C $genome > $RMout.TE.fa`;
`perl $cleanup_nested -in $RMout.TE.fa -threads $threads -minlen 80 -miniden 80 -cov 0.95 -blastplus $blast -iter 1 -maxcount $maxcount 2>/dev/null`;
for my $cat ("nested", "all", "redun") {
	`perl $count_nested -in $RMout.TE.fa.stat -cat $cat -mincov $mincov > $RMout.TE.fa.stat.$cat.sum`;
	for my $d (40, 30, 20, 10, 5) {
		`printf "\\n\\n" >> $RMout.TE.fa.stat.$cat.sum`;
		`perl $count_nested -in $RMout.TE.fa.stat -cat $cat -mincov $mincov -maxdiv $d >> $RMout.TE.fa.stat.$cat.sum`;
	}
}

$date=`date`;
chomp ($date);
print "$date\tEvaluation finished!\n\n";

