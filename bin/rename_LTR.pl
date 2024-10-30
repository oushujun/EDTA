#!/usr/bin/env perl
use warnings;
use strict;

#Usage: this script is developed to output gff3 format of intact LTRs and add classification info to LTR sequences that only have bare msu locus names
#perl rename_LTR.pl genome.fa target_sequence.fa LTR_retriever.defalse
#Shujun Ou (11/03/2019)

my $usage = "\n\tperl rename_LTR.pl genome.fa target_sequence.fa LTR_retriever.defalse\n\n";

my $genome = $ARGV[0];
my $seq = $ARGV[1]; #target seq
my $anno = $ARGV[2]; #annotation of the target seq
my $annotator = "EDTA";

open FA, "<$ARGV[0]" or die "ERROR: $usage";
open Seq, "<$seq" or die $usage;
open Anno, "<$anno" or die $usage;
open GFF, ">$seq.gff3" or die "ERROR: $!";

# read genome info and print the gff header
print GFF "##gff-version   3\n";
$/ = "\n>";
my $chr_info = '';
my %seq_flag;
my $j = 0;
while (<FA>){
	s/>//g;
	s/^\s+//;
	my ($id, $seq) = (split /\n/, $_, 2);
	$seq =~ s/\s+//g;
	my $len = length ($seq);
	print GFF "##sequence-region   $id 1 $len\n";
	$chr_info .= "#$id\n";
	$seq_flag{$id} = $j;
	$j++;
	}
print GFF "$chr_info";
$/ = "\n";

# read annotation info
my %anno;
my $i = "REPeaT";
while (<Anno>){
	next unless /motif/;
	my $nextline = <Anno>;

	my ($chr, $element_start, $element_end, $element_length, $lLTR_start, $lLTR_end, $lLTR_length, $rLTR_start, $rLTR_end, $rLTR_length, $seq_ID, $loc, $motif, $TSD, $lTSD, $rTSD, $IN, $sim, $strand, $supfam, $age);
	($loc, undef, $motif, $TSD, $lTSD, $rTSD, $IN, $sim, $strand, $supfam, undef, $age) = (split);
#chr03:19323647..19328532        pass    motif:TGCA      TSD:GTCGC       19323642..19323646      19328533..19328537      IN:19323854..19328325   99.03
#B73V4_ctg1:5489..14866  false   motif:CGGC      TSD:CTAT        5485..5488      14867..14870    IN:6808..13534  0.8923  +       Copia   LTR     4471309

	($lLTR_length, $rLTR_length) = ($1, $2) if $nextline =~ /lLTR: ([0-9a-z]+)\s+rLTR: ([0-9a-z]+)/i;
	next unless $lLTR_length ne "NA" and $rLTR_length ne "NA" and $lLTR_length > 0 and $rLTR_length > 0;
#        Adjust: NO      lLTR: NA        rLTR: NA
#        Adjust: 3' lLTR lLTR: 810       rLTR: 804
#        Adjust: NO      lLTR: 0 rLTR: 0

	($chr, $lLTR_start, $rLTR_end)=($1, $2, $3) if $loc=~/^(.*):([0-9]+)..([0-9]+)$/;
	($lLTR_start, $rLTR_end)=($rLTR_end, $lLTR_start) if $rLTR_end<$lLTR_start;
	($lLTR_end, $rLTR_start)=($1-1, $2+1) if $IN=~/^IN:([0-9]+)..([0-9]+)$/;
	$motif=~s/motif://gi;
	$TSD=~s/TSD://gi;
	$lTSD=~s/\.\./\t/;
	$rTSD=~s/\.\./\t/;
	$element_start=(split /\s+/, $lTSD)[0];
	$element_end=(split /\s+/, $rTSD)[1];
	my $id = "$chr:$lLTR_start..$rLTR_end#LTR/$supfam";
	if ($TSD eq "NA"){
		$element_start=$lLTR_start;
		$element_end=$rLTR_end;
		}
	my $chr_ori=$chr;
	my $gff = "$chr\t$annotator\trepeat_region\t$element_start\t$element_end\t.\t$strand\t.\tID=repeat_region$i\n";
	$gff .= "$chr\t$annotator\ttarget_site_duplication\t$lTSD\t.\t$strand\t.\tParent=repeat_region$i\n" unless $TSD eq "NA";
	$gff .= "$chr\t$annotator\tLTR/$supfam\t$lLTR_start\t$rLTR_end\t.\t$strand\t.\tID=$id;Parent=repeat_region$i;motif=$motif;tsd=$TSD;ltr_identity=$sim;seq_number=$seq_flag{$chr_ori}\n";
	$gff .= "$chr\t$annotator\tlong_terminal_repeat\t$lLTR_start\t$lLTR_end\t.\t$strand\t.\tParent=$id\n";
	$gff .= "$chr\t$annotator\tlong_terminal_repeat\t$rLTR_start\t$rLTR_end\t.\t$strand\t.\tParent=$id\n";
	$gff .= "$chr\t$annotator\ttarget_site_duplication\t$rTSD\t.\t$strand\t.\tParent=repeat_region$i\n" unless $TSD eq "NA";
	$gff .= "###\n";

	$anno{$loc} = [$supfam, $gff];
	}
close Anno;


# read target seq and print formatted names
$/ = "\n>";
my $k = 1;
while (<Seq>){
	s/>//g;
	my ($id, $seq) = (split /\n/, $_, 2);
	$id =~ s/\s+.*//;
	$id =~ s/#.*//;
	$seq =~ s/\s+//g;
	my ($chr, $lLTR_start, $rLTR_end);
	($chr, $lLTR_start, $rLTR_end) = ($1, $2, $3) if $id =~ /^(.*):([0-9]+)..([0-9]+)$/;
	($lLTR_start, $rLTR_end) = ($rLTR_end, $lLTR_start) if $rLTR_end < $lLTR_start;
	my ($anno, $gff) = ("unknown", "$chr\t$annotator\trepeat_region\t$lLTR_start\t$rLTR_end\t.\t?\t.\tID=repeat_region$k\n###\n");
	next unless defined $anno{$id}; #skip the seqs that are not in %anno (duplicated, heavily nested, or lack of LTR structure)
	if (defined $anno{$id}){
		($anno, $gff) = @{$anno{$id}};
		$gff =~ s/REPeaT/$k/g;
		}
	print ">$id#LTR/$anno\n$seq\n";
	print GFF "$gff";
	$k++;
	}
$/ = "\n";
close Seq;
close GFF;

