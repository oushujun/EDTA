#convert MGEScan-ltr candidate list to LTR_retriever input format
#Usage: perl convert.pl genome.ltrpos > genome.MGEScan.scn
#Author: Shujun OU (oushujun@msu.edu) 08/13/2016


#!/usr/bin/env perl -w
use strict;

print "#this summary file is created from MGEScan-ltr candidate list by convert_MGEScan.pl (Author: Shujun Ou, email: oushujun\@msu.edu 08/13/2016)
##start end len lLTR_str lLTR_end lLTR_len rLTR_str rLTR_end rLTR_len similarity seqid chr direction TSD lTSD rTSD motif superfamily family age(ya)\n";
my ($from, $to, $chr, $LTR_len, $lLTR_len, $rLTR_len, $lLTR_end, $rLTR_str, $similarity, $TSD, $motif, $direction, $seq_id);
while (<>){
	s/^\s+//;
	$from=$to=$LTR_len=$lLTR_len=$rLTR_len=$lLTR_end=$rLTR_str=$similarity=$direction="NA";
	($chr, $from, $lLTR_end, $rLTR_str, $to, $direction, $similarity)=(split)[0,2,3,4,5,6,10];
	#eg. Chr1    1       487160  487385  493666  493891  -       225     225     6281    98.7
	$from+=1; #1bp left-shift correction of start positions of the original output
	$lLTR_len=$lLTR_end-$from+1;
	$rLTR_len=$to-$rLTR_str+1;
	$LTR_len=$to-$from+1;
	$similarity="NA" unless $similarity ne '';
	$direction="NA" unless $direction ne '';
	print "$from $to $LTR_len $from $lLTR_end $lLTR_len $rLTR_str $to $rLTR_len $similarity $chr $chr $direction\n";
	}

