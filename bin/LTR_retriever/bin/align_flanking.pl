#!usr/bin/perl -w
use strict;

my $usage="
        Usage: perl align_sequence.pl aligned_length_cutoff(bp) similarity_cutoff(%) word_size Seq1.fa Seq2.fa
		perl align_flanking.pl a_cutoff s_cutoff w_size boundary_ctrl seq1.fa seq2.fa
        Example: perl align_sequence.pl 25 60 7 test1.fa test2.fa
        \n";
#Version
#	align_flanking.pl
#	align_flanking: To get the overall aligned length and identity of two sequences
#	Author: Shujun Ou (oushujun\@msu.edu), Department of Horticulture, Michigan State University, East Lansing, MI, 48824, USA
#	Version: 0.1 2014/01/27

my $a_cutoff=$ARGV[0]; #alignment length cutoff
my $s_cutoff=$ARGV[1]; #alignment identity cutoff
my $w_size=$ARGV[2]; #alignment word size
my $boundary_ctrl=$ARGV[3]; #boundary alignment control
my $blastplus='';
$blastplus=$ARGV[6] if defined $ARGV[6]; #path to the blast+ directory containing blastn

## Align the left and right boundaries
my ($bond, $seq3, $seq4, $left_align, $right_align, $left3, $left4, $right3, $right4, $boundary_aln);
$left_align=$right_align='None';
$boundary_aln="NA";
my $File3=$ARGV[4] or die "ERROR: $!";
my $File4=$ARGV[5] or die "ERROR: $!";
$File3=~s/^\s+//;
$File4=~s/^\s+//;
$seq3=(split /\\n/, $File3)[1];
$seq4=(split /\\n/, $File4)[1];
my $seq_l;
if (length($seq3)<=length($seq4)){
	$seq_l=length($seq3);
	}else{
	$seq_l=length($seq4);
	}
chomp ($seq3, $seq4);
$left3=substr $seq3, 0,10;
$left4=substr $seq4, 0,10;
$right3=substr $seq3, -10;
$right4=substr $seq4, -10;
$bond="l3:$left3 l4:$left4 r3:$right3 r4:$right4";
my @left3=split //, $left3;
my @left4=split //, $left4;
my @right3=split //, $right3;
my @right4=split //, $right4;
my ($i, $j, $a, $b, $align_cutoff);
$i=$j=$a=$b=0;
$align_cutoff=7;
foreach my $base (@left3){
	if ($base eq $left4[$i]){
		$a++;
		}
	$i++;
	}
foreach my $base (@right3){
	if ($base eq $right4[$j]){
		$b++;
		}
	$j++;
	}
if ($a>=$align_cutoff){
	$left_align="left";
	}
if ($b>=$align_cutoff){
	$right_align="right";
	}
$boundary_aln="$left_align,"."$right_align" if $boundary_ctrl==1;

## Align the flanking sequence
my @Blast=();
my $try=0;
while ($try<100){ #it's possible that sequence wrote in memory is rewritten by other programs and caused blast error, this step will try 100 times to guarantee the blast is run correctly
	@Blast=qx(bash -c '${blastplus}blastn -subject <(echo -e \"$File3\") -query <(echo -e \"$File4\") -evalue 1000 -word_size $w_size -dust no -outfmt 6' 2> /dev/null) if (defined $seq3 and defined $seq4);
#blast outfmt=6 looks like this:
#query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
#kasalath2       kasalath_chr01_91544_103018     88.00   25      3       0       7       31      1       25      5e-07   30.7
	last if $? == 0;
	$try++;
	}

my (%q_bank, %s_bank, $sim, $q_start, $q_end, $s_start, $s_end, $q_index, $s_index);
my ($length, $similarity, $m, $q_sim, $n, $s_sim, $s_e_align);
$sim=$q_start=$q_end=$s_start=$s_end=$similarity=$q_sim=$s_sim=$s_e_align=0.00;
$m=$n=$length=1;
$q_index=$s_index='';

unless (@Blast){ #if no blast result (not aligned), then goto the end of this script
	$length=0;
	goto End;
	}

foreach (@Blast){
	s/^\s+//;
	($sim, $q_start, $q_end, $s_start, $s_end)=(split)[2,6,7,8,9];
	$q_index="$q_index"."$q_start..$q_end;";
	$s_index="$s_index"."$s_start..$s_end;";
	if ($q_start>$q_end){
		($q_start, $q_end)=($q_end, $q_start);
		}
	if ($s_start>$s_end){
		($s_start, $s_end)=($s_end, $s_start);
		}
	if (10-$q_start>=4 && $s_end-(length($seq3)-10)>=4){ #15, 45
		$s_e_align=1; #head of query aligns with tail of subject (obsolete, original for A-D align)
		}

	for(my $i=$q_start; $i<=$q_end; $i++){
		unless (exists $q_bank{"q_$i"}){
			$q_bank{"q_$i"}=$sim;
			}
		}
	for (my $j=$s_start; $j<=$s_end; $j++){
		 unless (exists $s_bank{"s_$j"}){
			$s_bank{"s_$j"}=$sim;
			}
		}
	}

while ((my $key, my $value)=each (%q_bank)){
	$m++;
	$q_sim+=$value;
	}
while ((my $key, my $value)=each (%s_bank)){
	$n++;
	$s_sim+=$value;
	}

$length=($n+$m)/(2*$seq_l);
$similarity=($q_sim+$s_sim)/($m+$n);
$q_sim=sprintf("%.2f", $q_sim/$m);
$s_sim=sprintf("%.2f", $s_sim/$n);

End:
if ($length>=$a_cutoff && $similarity>=$s_cutoff){
	print "aligned";
	} else {
	print "not match";
	}
print "\tBoundary-align:$boundary_aln\t$bond\tHT-align:$s_e_align\tm:$m\tqsim:$q_sim\tn:$n\tssim:$s_sim\tqindex:$q_index\tsindex:$s_index";


