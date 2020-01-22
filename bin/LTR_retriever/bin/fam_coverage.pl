#!/usr/bin/env perl -w
use strict;

my $usage = "famcoverage.pl TE_lib RM_output genome_size_bp > TE_fam.size.list\n";

if (@ARGV < 3) {die "$usage\n";}

#This is for new RM output
#First version: Ning Jiang (jiangn@msu.edu)
#Last update: 03-21-2018 Shujun Ou (oushujun@msu.edu)
#Removed the dependency on seqstat from HMMER

my $whole = $ARGV[2]; 

open Lib, "<$ARGV[0]" || die "Can not open TE library: $ARGV[0]\n";

my (%RM_fam, %name, %rcnum, %subfam, %class, %subclass, %leng, %note, %cov, %ct, %full, %left, %right); 
$/ = "\n>";
while (<Lib>) {
	chomp;
	s/>//g;
	s/^\s+//;
	next if /^\s+$/;
	my ($name, $seq) = (split /\n/, $_, 2);
	$name = (split /\s+/, $name)[0];
	$seq =~ s/\s+//g;
	my $leng = length $seq;
	my ($fam, $subfam, $rcnum, $note, $class, $subclass);
	if ($name =~ /^((\S+)_(\S+))\#(\S+)\/(\S+)/) {
	    $fam = $1;
	    $rcnum = $2;
	    $note = $3;
	    $class = $4;
	    $subclass = $5;
	}
	elsif ($name =~ /^((\S+)_(\S+))\#(\S+)/) {
	    $fam = $1;
	    $rcnum = $2;
	    $note = $3;
	    $class = $4;
	    $subclass = "unknown";
	}
	elsif ($name =~ /^(\S+)\#(\S+)\/(\S+)/) {
	    $fam = $1;
	    $rcnum = $1;
	    $class = $2;
	    $subclass = $3;
	    $note = "NONE";
	}
	elsif ($name =~ /^(\S+)\#(\S+)/) {
	    $fam = $1;
	    $rcnum = $1;
	    if ($2 eq "MULETIR" || $2 eq "MULEtir"  ) {
	    $class = "DNA";
	    $subclass = "MULEtir";
	    $note = "NONE";
	    }
	    else {
		$class = $2;
	    $subclass = "unknown";
		$note = "NONE";
	    }
	}
	elsif  ($name =~ /^(\S+)/) {
	    $fam = $1;
	    $rcnum = $1;
	    $class = "unknown";
	    $subclass = "unknown";
	    $note = "NONE";
	}

 	if ($class ne "gene" && $class ne "RNA" && $class ne "plastid" && $class ne "vector" && $class ne "Bacterial") {
		my $RM_fam = $fam;
		$name{$RM_fam} = $name;
		$rcnum{$RM_fam} = $rcnum;
	        $subfam {$RM_fam} = '';
	        $class{$RM_fam} = $class;
	        $subclass{$RM_fam} = $subclass;
		$leng{$RM_fam} = $leng;
		$note{$RM_fam} = $note;
		$cov{$RM_fam} = 0;
		$ct{$RM_fam} = 0;
	        $full{$RM_fam} = 0;
		$left{$RM_fam} = 0;
		$right{$RM_fam} = 0;
		 }
	}
$/="\n";
close Lib;

open (MSK, "<$ARGV[1]") || die "Cannot open RM out $ARGV[1] file!\n";

my ($last_seq, $last_tail, $left_end, $right_end, $lt, $rt, $ft, $ctt, $ov, $covt) = ('', 0, 0, 0, 0, 0, 0, 0, 0, 0);
while (<MSK>) {
	s/^\s+//g;
	next unless /^[0-9]+/;
	s/[\(\)]//g;
	my ($seq, $head, $tail, $strand, $name, $subfam, $TE_head, $TE_tail, $TE_left)=(split)[4,5,6,8,9,10,11,12,13];
	$name =~ s/_INT-int$/_INT/;
	if (($seq ne $last_seq) or ($head > $last_tail)) {
		if ($strand eq "+") {
			$left_end = 1 if ($TE_head <= 20);
			$right_end = 1 if ($TE_left <= 20);
			}
		else {
			$right_end = 1 if ($TE_head <= 20);
			$left_end = 1 if ($TE_left <= 20);
			}
		$right_end = 0 if ($subfam =~ /tir/i || $name =~ /tir/i);
                }
	elsif ($head < $last_tail and $tail > $last_tail) {
		$ov = $head - $last_tail + 1; #the overlap between two overlapped annotations
		$head = $last_tail + 1;
		if ($tail - $head + 1 >= 20) {
			if ($strand eq "+") {
				$left_end = 1 if $TE_head + $ov <= 20;
				$right_end = 1 if ($TE_left <= 20);
				}
			else {
				$right_end = 1 if $TE_head + $ov <= 20;
				$left_end = 1 if ($TE_left <= 20);
				}
			$right_end = 0 if ($subfam =~ /tir/i || $name =~ /tir/i);
			}
		}

	if ($seq ne $last_seq){
		$cov{$name} += $tail - $head + 1; #cumulated length of the family
		$covt += $tail - $head + 1;
		$ct{$name} ++; #cumulated copy number of the family
		$ctt ++;
		}
	elsif ($head > $last_tail) {
		$cov{$name} += $tail - $head + 1;
		$covt += $tail - $head + 1;
		$ct{$name} ++;
		$ctt ++;
		}
	elsif ($tail > $last_tail) {
		$cov{$name} += $tail - $last_tail;
		$covt += $tail - $last_tail;
		$ct{$name} ++;
		$ctt ++;
	    }

	if ($left_end == 1 && $right_end == 1) {
		$full{$name} ++;
		$ft ++;
		}
	elsif ($left_end == 1) {
		$left{$name} ++;
		$lt ++;
		}
	elsif ($right_end == 1) {
		$right{$name} ++;
		$rt ++;
		}

	$last_tail = $tail;
	$last_seq = $seq;
	$left_end = 0;
	$right_end = 0;
	$ov = 0;
	}
close MSK;

printf "#RepeatMasker_entry TE_family Full_length Left_end_only Right_end_only Converted_copy_number Total_entries Total_length_in_bp Whole_genome_percentage Class Subclass Note\n";
foreach my $key (sort {$cov{$b} <=> $cov{$a}} (keys %cov)) {
	my $total = int ( $full{$key} + ($left{$key}+$right{$key})/2 + 0.5);
	my $coverage = sprintf("%.5f", $cov{$key}*100/$whole);
	print "$name{$key}\t$key\t$full{$key}\t$left{$key}\t$right{$key}\t$total\t$ct{$key}\t$cov{$key}\t$coverage\t$class{$key}\t$subclass{$key}\t$note{$key}\n";
}

my $total2 = int ($ft+($lt+$rt)/2 + 0.5);
my $coveage2 = sprintf("%.5f", $covt*100/$whole);
print "Summary: $ft\t$lt\t$rt\t$total2\t$ctt\t$covt\t$coveage2\n";

