#!/usr/bin/perl -w
use strict;
use FindBin;

my $usage = "\nFormat HelitronScanner fasta candidates with additional filterings
	perl format_helitronscanner_out.pl genome.fa
\n";

my $genome = $ARGV[0];
my $call_seq = "$FindBin::Bin/call_seq_by_list.pl";

my $ext_len = 10; #extend 10 bp on each end
my $tgt_ste_filter = 1; #1 will filter out candidate without AT or TT target site; 0 will not.
my $min_score = 12; #candidates with head and tail quality scores add up less than this will be discarded


die "HelitronScanner result files for the $genome is not found!\n$usage" unless -s $genome and -s "$genome.HelitronScanner.draw.rc.hel.fa";


open Hel, "cat $genome.HelitronScanner.draw.hel.fa $genome.HelitronScanner.draw.rc.hel.fa |" or die $usage;
open List, ">$genome.HelitronScanner.raw.ext.list" or die $usage;
my %hel;
while (<Hel>){
	next unless /^>/;
	chomp;
	s/>//;
	my ($loc, $dir, $len, $score) = (split)[0,1,2,4];
	my ($chr, $str, $end) = ($1, $2, $3) if $loc =~ /^(.*)_#SUB_([0-9]+)-([0-9]+)$/;

	$dir =~ s/[\[\]]+//g;
	#extend 5 bp on each end
	my ($new_str, $new_end);
	if ($dir eq "forward"){
		$new_str = $str - $ext_len;
		$new_end = $end + $ext_len;
		} else {
		$new_str = $str + $ext_len;
		$new_end = $end - $ext_len;
		}
	my $pos = "$chr:$new_str..$new_end";
	$score =~ s/scores=//;
	print List "$chr-$str-$end-$ext_len-$dir-$score\t$pos\n";
	}
close Hel;
close List;

`perl $call_seq $genome.HelitronScanner.raw.ext.list -C $genome > $genome.HelitronScanner.raw.ext.fa`;

open Hel2, "<$genome.HelitronScanner.raw.ext.fa" or die $usage;
open Out, ">$genome.HelitronScanner.filtered.tabout";
open Seq, ">$genome.HelitronScanner.filtered.fa";
print Out "#Chr\tStart\tEnd\tDirection\tLOC\tScore_head\tScore_tail\tTarget_site\t5'flank\t5'seq\t3'seq\t3'flank\n";
$/ = "\n>";
while (<Hel2>){
	chomp;
	s/>//g;
	my ($id, $seq) = (split /\n/, $_, 2);
	$seq =~ s/\s+//g;
	$id =~ s/^.*\|//;
	my ($chr, $str, $end, $ext_len, $dir, $score) = (split /-/, $id);
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
	next unless ($tgt_ste eq "TT" or $tgt_ste eq "AT") and $tgt_ste_filter == 1;
	next if $score_tot < $min_score;

	print Out "$chr\t$str\t$end\t$dir\t$loc\t$score_h\t$score_t\t$tgt_ste\t$flank5\t$seq5\t$seq3\t$flank3\n";
	print Seq ">$chr:$str..$end\t$dir|$tgt_ste|$score\n$helseq\n";
	}
close Hel2;
close Out;
close Seq;

