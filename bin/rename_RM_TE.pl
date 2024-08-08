#!/usr/bin/env perl
use warnings;
use strict;
# rename RepeatModeler consensus families, contributed by Zhigui Bao (#baozg)
# rename TEsorter families (Shujun Ou)
# Last modified: 04/19/2021

my $usage = "perl rename_RM_TE.pl RM_*/consensi.fa.classified > RepeatModeler.TE.raw.fa";
my $fasta = $ARGV[0];

my %hash=("DNA/hAT"=>"DNA/DTA",
          "TIR/hAT"=>"DNA/DTA",
          "DNA/CACTA"=>"DNA/DTC",
          "TIR/EnSpm_CACTA"=>"DNA/DTC",
          "DNA/PIF-Harbinger"=>"DNA/DTH",
          "TIR/PIF_Harbinger"=>"DNA/DTH",
          "DNA/Mutator"=>"DNA/DTM",
          "TIR/MuDR_Mutator"=>"DNA/DTM",
          "DNA/Tcl-Mariner"=>"DNA/DTT",
          "LTR" => "LTR/unknown",
          "LTR/Retrovirus" => "LTR/retrovirus",
          "LINE" => "LINE/unknown",
          "SINE" => "SINE/unknown",
          "SINE?" => "SINE/unknown",
          "Maverick" => "DNA/Maverick",
          "RC/Helitron"=>"DNA/Helitron",
          "Helitron"=>"DNA/Helitron",);


open FA, "<$fasta" or die "\nInput not found!\n$usage";
$/ = "\n>";
my $num = 0;
my %data;
while (<FA>){
	s/>//g;
	$num = sprintf("%08d", $num);
	my ($id, $seq) = (split /\n/, $_, 2);
	my $name = (split /\s/, $id)[0];
	$seq =~ s/\s+//g;
	my $class = $1 if $name =~ /^.*#(.*)/;
	next if $name =~ /mixture/i; #discard sequences that have mixture classifications
	#rename TE as unknown if $class info could not be retrieved
	my $tag = 0;
	foreach my $key (sort keys %hash){
		if ($class =~ /$key/i){
			print ">RM_${num}#$hash{$key}\n$seq\n";
			#print "$class\t$key\n";
			$tag = 1;
			last;
			}
		}
	if ($tag == 0){
		print ">RM_${num}#$class\n$seq\n";
		}
	$num +=1;
	}
close FA;
