#convert MGEScan 3.0.0 LTR candidates to LTR_retriever input format
#Usage: perl convert_MGEScan3.0.pl 
#convert.pl genome.ltrpos > genome.MGEScan.scn
#Author: Shujun Ou (shujun.ou.1@gmail.com) 04/19/2019


#!/usr/bin/env perl -w
use strict;

print "#this summary file is created from MGEScan 3.0.0 LTR candidate list by convert_MGEScan3.0.pl (Author: Shujun Ou, email: shujun.ou.1\@gmail.com 04/19/2019)
##start end len lLTR_str lLTR_end lLTR_len rLTR_str rLTR_end rLTR_len similarity seqid chr direction TSD lTSD rTSD motif superfamily family age(ya)\n";
while (<>){
	s/^$//;
	next if /^\s+$/;
	next if /^[0-9]+\-{5,}/;
	my ($chr, $from, $lLTR_end, $rLTR_str, $to, $direction, $similarity, $lLTR_len, $rLTR_len, $LTR_len, $lTSD, $rTSD, $lmotif1, $lmotif2, $rmotif1, $rmotif2);
	($chr, $from, $lLTR_end, $rLTR_str, $to, $direction, $lLTR_len, $rLTR_len, $LTR_len, $lTSD, $rTSD, $lmotif1, $lmotif2, $rmotif1, $rmotif2)=(split)[0,1,2,3,4,5,6,7,8,10,11,12,13,14,15];
	$chr=~s/\.fa_[0-9]+$//;
	my $TSD="$lTSD-$rTSD";
	my $motif = "${lmotif1}$lmotif2-${rmotif1}$rmotif2";
	$TSD="NA" if $TSD eq "---";
	$motif="NA" if $motif eq "-----";
	$similarity="NA";
	print "$from $to $LTR_len $from $lLTR_end $lLTR_len $rLTR_str $to $rLTR_len $similarity $chr $chr $direction $TSD $lTSD $rTSD $motif\n";
	}

