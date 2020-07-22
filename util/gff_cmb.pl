#!/usr/bin/env perl
use warnings;
use strict;

# Conbine gff files
# Shujun Ou (shujun.ou.1@gmail.com)

my $usage = "Usage: 
	cat *.gff3 | perl gff_cmb.pl - [structural|homology] > file_cmb.gff3 \n";
die "Please indicate if the specified GFF file is generated based on structural features or homology!\n$usage\n" unless defined $ARGV[1] and $ARGV[1] =~ /^structural$|^homology$/i;
my $method = $ARGV[1];
open GFF, "sort -sV -k5,5 $ARGV[0]|" or die $usage;

#defined sequence ontology
my %class = ("Cent/CentC" => "Cent", "Centro/tandem" => "Cent", "DNAauto/CACTA" => "TIR", "DNAauto/CACTG" => "TIR", "DNAauto/hAT" => "TIR", "DNAauto/Helitron" => "Helitron", "DNAauto/MLE" => "TIR", "DNAauto/MULE" => "TIR", "DNAauto/PILE" => "TIR", "DNAauto/POLE" => "TIR", "DNA/DTA" => "TIR", "DNA/DTC" => "TIR", "DNA/DTH" => "TIR", "DNA/DTM" => "TIR", "DNA/DTT" => "TIR", "DNA/Helitron" => "Helitron", "DNAnona/CACTA" => "TIR", "DNAnona/CACTG" => "TIR", "DNAnona/hAT" => "TIR", "DNAnona/Helitron" => "Helitron", "DNAnona/MLE" => "TIR", "DNAnona/MULE" => "TIR", "DNAnona/MULEtir" => "TIR", "DNAnona/PILE" => "TIR", "DNAnona/POLE" => "TIR", "DNAnona/Tourist" => "TIR", "DNAnona/unknown" => "DNATE", "Evirus/ERTBV-A" => "TE", "Evirus/ERTBV-B" => "TE", "Evirus/ERTBV-C" => "TE", "Evirus/ERTBV" => "TE", "knob/knob180" => "knob", "knob/TR-1" => "knob", "LINE/L1" => "LINE", "LINE/RTE" => "LINE", "LINE/unknown" => "LINE", "LTR/Copia" => "LTRRT", "LTR/CRM" => "LTRRT", "LTR/Gypsy" => "LTRRT", "LTR/Solo" => "LTRRT", "LTR/TRIM" => "LTRRT", "LTR/unknown" => "LTRRT", "MITE/DTA" => "MITE", "MITE/DTC" => "MITE", "MITE/DTH" => "MITE", "MITE/DTM" => "MITE", "MITE/DTT" => "MITE", "MITE/Stow" => "MITE", "MITE/Tourist" => "MITE", "Satellite/rice" => "satellite_DNA", "SINE/unknown" => "SINE", "subtelomere/4-12-1" => "subtelomere");

my %SO = (repeat_region => "SO:0000657", TE => "SO:0000101", DNATE => "SO:0000182", TIR => "SO:0000208", MITE => "SO:0000338", Helitron => "SO:0000544", retrotransposon => "SO:0000180", nonLTR => "SO:0000189", LINE => "SO:0000194", SINE => "SO:0000206", LTRRT => "SO:0000186", target_site_duplication => "SO:0000434", primer_binding_site => "SO:0005850", long_terminal_repeat => "SO:0000286", Low_complexity => "SO:0001005", "rDNA/spacer" => "SO:0001860", telomeric_repeat => "SO:0001496", subtelomere => "SO:0001997", Cent => "SO:0001797", satellite_DNA => "SO:0000005");

#print out gff3 head
my $date=`date -u`;
chomp ($date);
print "##gff-version 3\n##date $date
##seqid source repeat_class/superfamily start end sw_score strand phase attributes\n";

my $i = 0; #annotation ID
while (<GFF>){
	chmod;
	next if /^#/;
	my ($SW_score, $div, $iden, $chr, $chr_len, $element_start, $element_end, $element_length, $left_len, $strand, $TE_ID, $TE_class);
	(
	($SW_score, $div, $chr, $element_start, $element_end, $left_len, $strand, $TE_ID, $TE_class)=(split)[0,1,4,5,6,7,8,9,10];

	my ($chr, $TE_class, $element_start, $element_end, $score, $dir, $info) = (split)[0,2,3,4,5,6,8];
	my $class = "undef";
	$class = "LTR" if $type =~ /LTR/i or $type =~ /long_terminal_repeat/i or $type =~ /target_site_duplication/i;
	$class = "TIR" if $type =~ /DT/ or ($type =~ /DNA|MITE/ and $type !~ /Helitron|DHH/);
	$class = "Cent" if $type =~ /Cent/i;
	$class = "knob" if $type =~ /knob/i;
	$class = "LINE" if $type =~ /LINE|RIL/i;
	$class = "SINE" if $type =~ /SINE|RIS/i;
	$class = "rDNA" if $type =~ /rDNA/i;
	$class = "SAT" if $type =~ /SAT/i;
	$class = "subtelomere" if $type =~ /subtelomere/i;
	$class = "Helitron" if $type =~ /Helitron|DHH/i;
	$class = $1 if $type =~ /^(.*)\/.*/ and $1 !~ /DNA|MITE/i;
	$info =~ s/ID=TE_annot_[0-9]+;//; #rename annotation id based on input order
	$info =~ s/;Method=.*//; #remove method info at the end of the annotation
	next if $type eq "repeat_region";
	next unless defined $id and defined $start;
	print "$id\t$start\t$end\t$class\t$type\t$score\t$dir\tID=TE_annot_$i;$info;Method=$method\n";
	$i++;
	}

