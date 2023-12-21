#!/usr/bin/env perl
use warnings;
use strict;
use FindBin;

my $usage = "\nFilter HelitronScanner fasta candidates
	perl format_helitronscanner_out.pl -genome genome.fa [options]
		-sitefilter	[0|1]	1 will filter out candidate without AT or TT target site (default); 0 will not.
		-minscore	[int]	Candidates with head and tail quality scores add up less than this will be discarded. Default: 12
		-keepshorter	[0|1]	1 will keep the shorter possible when multi 5' end presents (default); 0 will not; 2 will keep all possible 5' ends.
		-extlen	[int]	Length of flanking sequence for blast and output. Default: 30 (bp)
		-extout	[0|1]	Output original sequence (0, default) or extended (1) sequence.
		-h|-help	Display this help messege and exit.
\n";

my $genome = '';
my $ext_len = 30; #extend 30 bp on each end
my $ext_out = 0; #Output original sequence (0, default) or extended (1) sequence.
my $tgt_ste_filter = 1; #1 will filter out candidate without AT or TT target site; 0 will not.
my $min_score = 12; #candidates with head and tail quality scores add up less than this will be discarded
my $keep_shorter = 1; #0 will keep the default 5' end; 1 will keep the shorter possible when multi 5' end presents (default); 2 will output all 5' ends.
my $call_seq = "$FindBin::Bin/call_seq_by_list.pl";

my $k=0;
foreach (@ARGV){
        $genome = $ARGV[$k+1] if /^-genome$/i;
	$ext_len = $ARGV[$k+1] if /^-extlen$/i;
	$ext_out = $ARGV[$k+1] if /^-extout$/i;
	$tgt_ste_filter = $ARGV[$k+1] if /^-sitefilter$/i;
	$min_score = $ARGV[$k+1] if /^-minscore$/i;
	$keep_shorter = $ARGV[$k+1] if /^-keepshorter$/i;
	die $usage if /^-h$|^-help$/i;
	$k++;
	}

die "HelitronScanner result files for the $genome is not found!\n$usage" unless -e $genome and -e "$genome.HelitronScanner.draw.rc.hel.fa";


open Hel, "cat $genome.HelitronScanner.draw.hel.fa $genome.HelitronScanner.draw.rc.hel.fa |" or die $usage;
open List, ">$genome.HelitronScanner.raw.ext.list" or die $usage;
my %hel;
while (<Hel>){
	next unless /^>/;
	chomp;
	s/>//;
	my ($loc, $dir, $len, $score, $alt5) = (split /\s+/, $_, 6)[0,1,2,4,5];
	my ($chr, $str, $end) = ($1, $2, $3) if $loc =~ /^(.*)_#SUB_([0-9]+)-([0-9]+)$/;

	$dir =~ s/[\[\]]+//g;
	$alt5 =~ s/Multi_5'_ends://;
	$score =~ s/scores=//;

	# store default 5' end
	my %line;
	$line{"$str-$end-$score"} = '' if $keep_shorter != 1;

	# store alternative 5' ends
	if ($alt5 ne ''){
		my $short5 = $str;
		while ($alt5 =~ s/([0-9]+):([0-9]+)//){
			my ($test5, $test_score) = ($1, $2);
			$score =~ s/^([0-9]+):/$test_score:/;
			$line{"$test5-$end-$score"} = '' if $keep_shorter == 2; #keep all 5' ends
			$short5 = $test5 if abs($short5 - $end) > abs($test5 - $end);
			}
		$line{"$short5-$end-$score"} = '' if $keep_shorter == 1; #keep shorter 5' ends
		} 
	elsif ($keep_shorter == 1){
		$line{"$str-$end-$score"} = '';
		}

	# process all 5' ends
	foreach (keys %line){
		my ($this_str, $this_end, $this_score) = (split /\-/, $_);
		#extend $ext_len bp on each end
		my ($new_str, $new_end);
		if ($dir eq "forward"){
			$new_str = $this_str - $ext_len;
			$new_end = $this_end + $ext_len;
			} else {
			$new_str = $this_str + $ext_len;
			$new_end = $this_end - $ext_len;
			}
		my $pos = "$chr:$new_str..$new_end";
		print List "$chr-$this_str-$this_end-$ext_len-$dir-$this_score\t$pos\n";
		}
	}
close Hel;
close List;

## Get extended fasta seq
`perl $call_seq $genome.HelitronScanner.raw.ext.list -C $genome > $genome.HelitronScanner.raw.ext.fa`;

open Hel2, "<$genome.HelitronScanner.raw.ext.fa" or die $usage;
open Out, ">$genome.HelitronScanner.filtered.tabout";
if ($ext_out eq 1){
	open Seq, ">$genome.HelitronScanner.filtered.ext.fa";
	} else {
	open Seq, ">$genome.HelitronScanner.filtered.fa";
	}
print Out "#Chr\tStart\tEnd\tDirection\tLOC\tScore_head\tScore_tail\tTarget_site\t5'flank\t5'seq\t3'seq\t3'flank\n";
$/ = "\n>";
while (<Hel2>){
	chomp;
	s/>//g;
	my ($id, $seq) = (split /\n/, $_, 2);
	$seq =~ s/\s+//g;
	$id =~ s/^.*\|//;
	my ($chr, $str, $end, $ext_len, $dir, $score);
	($chr, $str, $end, $ext_len, $dir, $score) = ($1, $2, $3, $4, $5, $6) if $id =~ /(.*)\-([0-9]+)\-([0-9]+)\-([0-9]+)\-(forward|reverse)\-([0-9:]+)/;
	next unless defined $chr;
	my $loc = "$chr:$str..$end";
	if ($dir eq "forward"){
		$dir = "+";
		}
	elsif ($dir eq "reverse") {
		$dir = "-";
		}
	my ($score_h, $score_t, $score_tot) = (0, 0, 0);
	($score_h, $score_t) = (split /:/, $score);
	$score_tot = $score_h + $score_t;

	my ($flank5, $flank3, $seq5, $seq3, $tgt_ste, $helseq) = ('','','','','','');
	$flank5 = substr $seq, 0, $ext_len;
	$flank3 = substr $seq, -$ext_len;
	$seq5 = substr $seq, $ext_len, 30;
	$seq3 = substr $seq, -($ext_len+30), 30;
	$tgt_ste = (substr $flank5, -1).(substr $flank3, 0, 1);
	$helseq = substr $seq, $ext_len, -$ext_len;

	#filter out candidates based on target site and score
	$tgt_ste = uc $tgt_ste;
	if ($tgt_ste_filter == 1) {
		next unless ($tgt_ste eq "TT" or $tgt_ste eq "AT");
		}
	next if $score_tot < $min_score;

	print Out "$chr\t$str\t$end\t$dir\t$loc\t$score_h\t$score_t\t$tgt_ste\t$flank5\t$seq5\t$seq3\t$flank3\n";
	print Seq ">$chr:$str..$end\t$dir|$tgt_ste|$score|$ext_len\n$helseq\n" if $ext_out eq 0;
	print Seq ">$chr:$str..$end\t$dir|$tgt_ste|$score|$ext_len\n$seq\n" if $ext_out eq 1;
	}
close Hel2;
close Out;
close Seq;
