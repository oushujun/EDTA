#! /usr/bin/perl

$usage = "fam_summary.pl TE_fam.size.list genome_size_bp > TE_fam.sum.txt\n";

if (@ARGV < 2) {die "$usage\n";}

#This is for new RM output

$whole = $ARGV[1]; 

open(MSK, "sort -k 10,10 -k 11,11 -k 12,12 $ARGV[0] |") || die "cannot open $ARGV[0]";


print " full    left    right   total    all      bps     percent     class      subclass        note\n";
while (<MSK>) {
	if (/^\s*\S+\s+\S+\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)\s*/) {
	    $class = $7;
	    $subc = $8;
	    $note = $9;
	    if ($class ne $lclass && $bps > 0) {
		printf "%06d  %06d  %06d  %06d  %07d  %09d  %.5f\t%-10s  %-15s  %s\n", $full,$left,$right,$total,$all,$bps,($bps*100/$whole),$lclass,$lsubc,$lnote;
		$full = $1;
		$left = $2;
		$right = $3;
		$total = $4;
		$all = $5;
		$bps = $6;
	    }
	    elsif ($subc ne $lsubc && $bps > 0) {
		printf "%06d  %06d  %06d  %06d  %07d  %09d  %.5f\t%-10s  %-15s  %s\n", $full,$left,$right,$total,$all,$bps,($bps*100/$whole),$lclass,$lsubc,$lnote;
		$full = $1;
		$left = $2;
		$right = $3;
		$total = $4;
		$all = $5;
		$bps = $6;
	    }
	    else {
		$full = $1 + $full;
		$left = $2 + $left;
		$right = $3 + $right;
		$total = $4 + $total;
		$all = $5 + $all;
		$bps = $6 + $bps;
	    }
	    $lclass = $class;
	    $lsubc = $subc;
	    $lnote = $note;
	}

}
close MSK;

printf "%06d  %06d  %06d  %06d  %07d  %09d  %.5f\t%-10s  %-15s  %s\n", $full,$left,$right,$total,$all,$bps,($bps*100/$whole),$lclass,$lsubc,$lnote;
