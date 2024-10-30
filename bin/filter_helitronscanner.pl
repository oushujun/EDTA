#!/usr/bin/env perl
use warnings;
use strict;
use FindBin;
use threads;
use Thread::Queue;

my $usage = "\nFilter HelitronScanner fasta candidates
	perl filter_helitronscanner_out.pl -genome genome.fa [options]
		-sitefilter	[0|1]	1 will filter out candidate without AT or TT target site (default); 0 will not.
		-minscore	[int]	Candidates with head and tail quality scores add up less than this will be discarded. Default: 12
		-keepshorter	[0|1]	1 will keep the shorter possible when multi 5' end presents (default); 0 will not.
		-extlen	[int]	Length of flanking sequence for blast and output. Default: 30 (bp)
		-miniden	[int]	Minimum identity for flanking sequence alignment. Default: 80 (%)
		-mincov	[float]	Minimum coverage for flanking sequence alignment that counts as full match. Default: 0.8
		-blastplus	[path]	Path to the blastn program. Defalut: read from \$ENV
		-t|-threads	[int]	Number of threads to run this program. Default: 4
		-h|-help	Display this help messege and exit.
\n";

my $genome = '';
my $ext_len = 30; #extend 30 bp on each end
my $tgt_ste_filter = 1; #1 will filter out candidate without AT or TT target site; 0 will not.
my $min_score = 12; #candidates with head and tail quality scores add up less than this will be discarded
my $keep_shorter = 1; #1 will keep the shorter possible when multi 5' end presents (default); 0 will not.
my $min_iden = 80; #minimum identity for flanking sequence alignment (%)
my $min_cov = 0.8; #minimum coverage for flanking sequence alignment that counts as full match
my $blastplus = ''; #path to the blastn program
my $threads = 4; #threads to run this program
my $call_seq = "$FindBin::Bin/call_seq_by_list.pl";

my $k=0;
foreach (@ARGV){
        $genome = $ARGV[$k+1] if /^-genome$/i;
	$ext_len = $ARGV[$k+1] if /^-extlen$/i;
	$tgt_ste_filter = $ARGV[$k+1] if /^-sitefilter$/i;
	$min_score = $ARGV[$k+1] if /^-minscore$/i;
	$keep_shorter = $ARGV[$k+1] if /^-keepshorter$/i;
	$min_iden = $ARGV[$k+1] if /^-miniden$/i;
	$min_cov = $ARGV[$k+1] if /^-mincov$/i;
	$blastplus=$ARGV[$k+1] if /^-blastplus$/i;
	$threads=$ARGV[$k+1] if /^-t$|^-threads$/i;
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

	# get shorter coordinates for candidates with alternative 5'end available.
	if ($keep_shorter == 1 and $alt5 ne ''){
		my $short5 = $str;
		while ($alt5 =~ s/([0-9]+):[0-9]+//){
			my $test5 = $1;
			$short5 = $test5 if abs($short5 - $end) > abs($test5 - $end);
			}
		$str = $short5;
		}

	#extend $ext_len bp on each end
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

## Get extended fasta seq
`perl $call_seq $genome.HelitronScanner.raw.ext.list -C $genome > $genome.HelitronScanner.raw.ext.fa`;

open Hel2, "<$genome.HelitronScanner.raw.ext.fa" or die $usage;
open Out, ">$genome.HelitronScanner.filtered.cov${min_cov}iden${min_iden}.tabout";
open Seq, ">$genome.HelitronScanner.filtered.fa";
#open Seq, ">$genome.HelitronScanner.filtered.cov${min_cov}iden${min_iden}.fa";
print Out "#Decision\t5'count\t3'count\tChr\tStart\tEnd\tDirection\tLOC\tScore_head\tScore_tail\tTarget_site\t5'flank\t5'seq\t3'seq\t3'flank\n";

## Store sequence information
my @FA;
$/ = "\n>";
while (<Hel2>){
	chomp;
	s/>//g;
	my ($id, $seq) = (split /\n/, $_, 2);
	$id =~ s/^.*\|//;
	$seq =~ s/\s+//g;
	$seq = uc $seq;
	push @FA, [$id, $seq];
	}
$/ = "\n";
close Hel2;

## multi-threading using queue, put candidate LTRs into queue for parallel computation
my $queue=Thread::Queue->new();
my $i=0;
while ($i<=$#FA) {
	last unless defined $FA[$i]->[0];
	my ($name, $seq)=@{$FA[$i]}[0,1];
	my ($chr, $seq_start, $seq_end, $ltr_start, $ltr_end);
	$queue->enqueue([$name, $seq]);
	$i++;
	}

## initiate a number of worker threads
my @threads=();
foreach (1..$threads){
	push @threads,threads->create(\&filter);
	}
foreach (@threads){
	$queue->enqueue(undef);
	}
foreach (@threads){
	$_->join();
	}
close Out;
close Seq;

## fixing the formatting error created by simutaniously writing the same file
`perl -i -nle 's/>/\\n>/g unless /^>/; print \$_' $genome.HelitronScanner.filtered.fa`;




## subrotine for helitron candidate analyses
sub filter(){
	while (defined($_ = $queue->dequeue())){
	my ($id, $seq)=(@{$_}[0], @{$_}[1]);
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

	# filter out candidates based on target site and score
	$tgt_ste = uc $tgt_ste;
	next unless ($tgt_ste eq "TT" or $tgt_ste eq "AT") and $tgt_ste_filter == 1;
	next if $score_tot < $min_score;

	# filter out candidates if flanking sequences are simple repeat
	my ($ssr_flank5, $ssr_flank3) = ('NA', 'NA');
	($ssr_flank5, $ssr_flank3) = (&count_base($flank5), &count_base($flank3));
	next if $ssr_flank5 eq 'true' or $ssr_flank3 eq 'true';

	# filter out candidates based on repetitiveness of flanking sequence
	my $decision="true";

	# count copy number of the 5' end
	my $end5_repeat = "false";
	my $end5 = ">end5\\n$flank5"."$seq5";
	my $end5_len = length "$flank5"."$seq5";
	my $exec="${blastplus}blastn -subject $genome -query <(echo -e \"$end5\") -outfmt 6 -word_size 7 -evalue 1e-5 -dust no";
	my @blast_end5=();
	my $try=0;
	while ($try<10){ #try 10 times to guarantee the blast is run correctly
		@blast_end5=qx(bash -c '$exec' 2> /dev/null) if defined $end5;
		last if $? == 0;
		$try++;
		}
	my $end5_count=0;
	foreach (@blast_end5){
		my ($iden, $len) = (split)[2,3];
		$end5_count++ if $iden >= $min_iden and $len >= $end5_len * $min_cov;
		($end5_repeat = "true", $decision = "false") if $end5_count > 1;
		}
	$end5_count = 'NA' unless @blast_end5 > 0;

	# count copy number of the 3' end
	my $end3_repeat = "false";
	my $end3 = ">end3\\n$seq3"."$flank3";
	my $end3_len = length "$flank3"."$seq3";
	$exec="${blastplus}blastn -subject $genome -query <(echo -e \"$end3\") -outfmt 6 -word_size 7 -evalue 1e-5 -dust no";
	my @blast_end3=();
	$try=0;
	while ($try<10){
		@blast_end3=qx(bash -c '$exec' 2> /dev/null) if defined $end3;
		last if $? == 0;
		$try++;
		}
#print "$id\t$end3_len\n@blast_end3\n";
	my $end3_count=0;
	if ($#blast_end3>0){
		foreach (@blast_end3){
			my ($iden, $len) = (split)[2,3];
			$end3_count++ if $iden >= $min_iden and $len >= $end3_len * $min_cov;
			($end3_repeat = "true", $decision = "false") if $end3_count > 1;
			}
		}
	$end3_count = 'NA' unless @blast_end3 > 0;

	# count copy number of the 5' and 3' flanking. If $count>=1, then this candidate locates at a TE and should be a true helitron
	if ($end5_repeat eq "true" and $end3_repeat eq "true"){
		my $flank = ">flank\\n$flank5"."$flank3";
		my $flank_len = length "$flank5"."$flank3";
		$exec="${blastplus}blastn -subject $genome -query <(echo -e \"$flank\") -outfmt 6 -word_size 7 -evalue 1e-5 -dust no";
		my @blast_flank=();
		$try=0;
		while ($try<10){
			@blast_flank=qx(bash -c '$exec' 2> /dev/null) if defined $flank;
			last if $? == 0;
			$try++;
			}
		if ($#blast_flank>0){
			my $count=0;
			foreach (@blast_flank){
				my ($iden, $len) = (split)[2,3];
				$count++ if $iden >= $min_iden and $len >= $flank_len * $min_cov;
				}
			$decision = "true" if $count >= 1;
			}
		}

	# if flanking blast failed, mark the candidate as false
	$decision = "false" if $end5_count eq 'NA' or $end3_count eq 'NA';

	print Out "$decision\t$end5_count\t$end3_count\t$chr\t$str\t$end\t$dir\t$loc\t$score_h\t$score_t\t$tgt_ste\t$flank5\t$seq5\t$seq3\t$flank3\n";
	print Seq ">$chr:$str..$end\t$dir|$tgt_ste|$score\n$helseq\n" if $decision eq "true";
	}
	}

# determine if the given sequence is simple repeat
sub count_base(){
	my $seq = lc $_[0];
	$seq =~ s/[nx]+//gi;
	my @base = ('a', 't', 'c', 'g');
	my @count = map { $_ = () = ($seq =~ /$_/gi) } @base; #base count
	@count = (sort { $b<=>$a } @count); #reverse sort
	my $dominant_base = ($count[0] + $count[1])/length $seq;
	my $repeat = "false";
	$repeat = "true" if $dominant_base >= 0.85;
	return $repeat;
	}
