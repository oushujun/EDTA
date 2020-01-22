##Description: Identify TE type (LTR or not LTR), LTR family (copia or gypsy), and strand of the input hmmsearch table
##		hmm profile categories were made based on the rice TE library "rice6.8.liban"
##Usage:       perl annotate_TE.pl hmmsearch.tbl > TE.family-strand.tbl
##		To obtain the hmmsearch.tbl, the recommend hmmsearch command is:
#		$ hmmsearch --tblout hmmsearch.tbl --notextw --cpu 5 -E 0.05 --domE 0.05 --noali TE.hmm input.aa.fa > hmmsearch.scn
##Author:      Shujun Ou (oushujun@msu.edu) 
##Versions:    v1.0	06/16/2016


#!/usr/bin/env perl -w
use strict;

my %notLTR=qw(PF10551.6 notLTR PF13966.3 notLTR PF04434.14 notLTR PF03108.12 notLTR PF13359.3 notLTR PF00872.15 notLTR PF05699.11 notLTR PF13963.3 notLTR PF04827.11 notLTR PF05970.11 notLTR PF03101.12 notLTR PF13960.3 notLTR PF02992.11 notLTR PF14214.3 notLTR PF12776.4 notLTR PF13952.3 notLTR PF03372.20 notLTR PF14372.3 notLTR PF14303.3 notLTR PF13604.3 notLTR PF01498.15 notLTR PF01609.18 notLTR PF13837.3 notLTR PF13538.3 notLTR PF02892.12 notLTR PF13384.3 notLTR PF06839.9 notLTR PF03004.11 notLTR PF04937.12 notLTR PF13191.3 notLTR PF02902.16 notLTR PF12728.4 notLTR PF14291.3 notLTR PF14529.3 notLTR PF10446.6 notLTR PF14111.3 notLTR PF13873.3 notLTR PF13245.3 notLTR PF04147.9 notLTR PF13542.3 notLTR PF05285.9 notLTR PF00376.20 notLTR PF00249.28 notLTR PF00646.30 notLTR PF12116.5 notLTR PF09322.7 notLTR PF00437.17 notLTR PF13401.3 notLTR PF06869.9 notLTR PF13857.3 notLTR PF04967.9 notLTR PF00564.21 notLTR PF13412.3 notLTR PF12796.4 notLTR PF06072.8 notLTR);

my %LTR=qw(PF00665.23 LTR PF03732.14 LTR PF14223.3 LTR PF07727.11 LTR PF13650.3 LTR PF13975.3 LTR PF13976.3 LTR PF08284.8 LTR PF00385.21 LTR PF13683.3 LTR PF09337.7 LTR PF13961.3 LTR PF09668.7 LTR PF14244.3 LTR PF13917.3 LTR PF04195.9 LTR PF03578.12 LTR PF14787.3 LTR PF02160.12 LTR PF05754.11 LTR PF04094.11 LTR PF16025.2 LTR PF11835.5 LTR PF12435.5 LTR PF05585.9 LTR PF05377.8 LTR PF04582.9 LTR PF12353.5 LTR PF06882.9 LTR PF08847.8 LTR PF10046.6 LTR PF17150.1 LTR PF12384.5 LTR PF01519.13 LTR PF00077.17 LTR PF15450.3 LTR PF11180.5 LTR PF08149.8 LTR PF07889.9 LTR PF12425.5 LTR PF07197.9 LTR PF05837.9 LTR PF04513.9 LTR PF03154.12 LTR PF04111.9 LTR PF04102.9 LTR PF06827.11 LTR);

my %gypsy=qw(PF13650.3 gypsy PF13975.3 gypsy PF08284.8 gypsy PF00385.21 gypsy PF09337.7 gypsy PF13456.3 gypsy PF04195.9 gypsy PF03578.12 gypsy PF00075.21 gypsy PF05754.11 gypsy PF04094.11 gypsy PF05585.9 gypsy PF06882.9 gypsy PF00077.17 gypsy PF10536.6 gypsy PF00078.24 gypsy PF09668.7 gypsy PF02160.12 gypsy PF16025.2 gypsy PF11835.5 gypsy PF04582.9 gypsy PF05377.8 gypsy PF17150.1 gypsy PF10046.6 gypsy PF12384.5 gypsy PF12353.5 gypsy PF07197.9 gypsy PF15450.3 gypsy PF04513.9 gypsy PF08149.8 gypsy PF03732.14 gypsy PF00665.23 gypsy);

my %copia=qw(PF07727.11 copia PF13961.3 copia PF14244.3 copia PF03564.12 copia PF13976.3 copia PF14223.3 copia PF08847.8 copia PF14787.3 copia PF13917.3 copia);

my %TE;
my %hmm;
open File, "<$ARGV[0]" or die "ERROR: $ARGV[0]. $!";
while (<File>){
	next if /^#/;
	s/^\s+//;
	my ($TE, $hmm)=(split)[0,3];
	my ($decision, $family, $strand)=("NA", "unknown", "?");

#Identify strand of the current entry
	$TE=~s/^(.*)\|(.*aa[1-3])$/$1/;
	$strand=$2;
	if ($strand=~/rev/i){
		$strand="-";
		} else {
		$strand="+";
		}

#judge the current entry
	$decision="LTR" if exists $LTR{$hmm};
	$decision="notLTR" if exists $notLTR{$hmm};

#idenfity family of the current entry
	$family="Gypsy" if exists $gypsy{$hmm};
	$family="Copia" if exists $copia{$hmm};

#create the initial TE profile
#	Decision	gypsyhmm#	copiahmm#	+strand#	-strand#	hmmmatch
	$TE{$TE}=[$decision, "0", "0", "0", "0", $hmm] unless exists $TE{$TE};

#update judge of the TE
	if (defined $TE{$TE}[0]){
		$TE{$TE}[5].=" $hmm";
		$TE{$TE}[0]="mixture" if $TE{$TE}[0] ne $decision;
		}

#count hmm type of the TE
	$TE{$TE}[1]++ if $family eq "Gypsy";
	$TE{$TE}[2]++ if $family eq "Copia";

#count strand type of the TE
	$TE{$TE}[3]++ if  $strand eq "+";
	$TE{$TE}[4]++ if  $strand eq "-";
	}
close File;

print "#TE\tSuperfamily\tFamily\tStrand\thmmmatchs\n";
foreach (keys %TE){
#Final decision of family and strand based on occurrance frequency
	my ($gypsy, $copia,$pos_strand, $neg_strand)=@{$TE{$_}}[1,2,3,4];
	my ($family, $strand)=("unknown", "?");
	$family="Gypsy" if ($gypsy>$copia and $copia/$gypsy<0.3);
	$family="Copia" if ($copia>$gypsy and $gypsy/$copia<0.9); #improvement suggested by @zhangrengang
	$strand="+" if ($pos_strand>$neg_strand and $neg_strand/$pos_strand<0.4);
	$strand="-" if ($neg_strand>$pos_strand and $pos_strand/$neg_strand<0.4);

#TE	Decision	family	strand	hmmmatch
	print "$_\t$TE{$_}[0]\t$family\t$strand\t$TE{$_}[5]\n";
	}

