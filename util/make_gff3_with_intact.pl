#!/usr/bin/env perl
use warnings;
use strict;
#Shujun Ou (shujun.ou.1@gmail.com) 11/15/2019

my $usage = "\n\tperl make_gff_with_intact.pl EDTA.intact.fa\n\n";
die $usage unless -s $ARGV[0];
my $date=`date -u`;
chomp ($date);
open ID, "grep \\> $ARGV[0] | sort -sV -k1,1 |" or die "ERROR: $!";
open GFF, ">$ARGV[0].gff3" or die "ERROR: $!";
print GFF "##gff-version 3\n##date $date
##seqid source repeat_class/superfamily start end score strand phase attributes\n";
##Chromosome Annotator Repeat_class/superfamily Start End Score Strand Phase Repeat_famliy\n";

#defined sequence ontology
my %class = ("Cent/CentC" => "Cent", "Centro/tandem" => "Cent", "DNAauto/CACTA" => "TIR", "DNAauto/CACTG" => "TIR", "DNAauto/hAT" => "TIR", "DNAauto/Helitron" => "Helitron", "DNAauto/MLE" => "TIR", "DNAauto/MULE" => "TIR", "DNAauto/PILE" => "TIR", "DNAauto/POLE" => "TIR", "DNA/DTA" => "TIR", "DNA/DTC" => "TIR", "DNA/DTH" => "TIR", "DNA/DTM" => "TIR", "DNA/DTT" => "TIR", "DNA/Helitron" => "Helitron", "DNAnona/CACTA" => "TIR", "DNAnona/CACTG" => "TIR", "DNAnona/hAT" => "TIR", "DNAnona/Helitron" => "Helitron", "DNAnona/MLE" => "TIR", "DNAnona/MULE" => "TIR", "DNAnona/MULEtir" => "TIR", "DNAnona/PILE" => "TIR", "DNAnona/POLE" => "TIR", "DNAnona/Tourist" => "TIR", "DNAnona/unknown" => "DNATE", "Evirus/ERTBV-A" => "TE", "Evirus/ERTBV-B" => "TE", "Evirus/ERTBV-C" => "TE", "Evirus/ERTBV" => "TE", "knob/knob180" => "knob", "knob/TR-1" => "knob", "LINE/L1" => "LINE", "LINE/RTE" => "LINE", "LINE/unknown" => "LINE", "LTR/Copia" => "LTRRT", "LTR/CRM" => "LTRRT", "LTR/Gypsy" => "LTRRT", "LTR/Solo" => "LTRRT", "LTR/TRIM" => "LTRRT", "LTR/unknown" => "LTRRT", "MITE/DTA" => "MITE", "MITE/DTC" => "MITE", "MITE/DTH" => "MITE", "MITE/DTM" => "MITE", "MITE/DTT" => "MITE", "MITE/Stow" => "MITE", "MITE/Tourist" => "MITE", "Satellite/rice" => "satellite_DNA", "SINE/unknown" => "SINE", "subtelomere/4-12-1" => "subtelomere");

my %SO = (repeat_region => "SO:0000657", TE => "SO:0000101", DNATE => "SO:0000182", TIR => "SO:0000208", MITE => "SO:0000338", Helitron => "SO:0000544", retrotransposon => "SO:0000180", nonLTR => "SO:0000189", LINE => "SO:0000194", SINE => "SO:0000206", LTRRT => "SO:0000186", target_site_duplication => "SO:0000434", primer_binding_site => "SO:0005850", long_terminal_repeat => "SO:0000286", Low_complexity => "SO:0001005", "rDNA/spacer" => "SO:0001860", telomeric_repeat => "SO:0001496", subtelomere => "SO:0001997", Cent => "SO:0001797", satellite_DNA => "SO:0000005");


my $annotator="EDTA";
my $i = 0; #annotation count
while (<ID>){
	s/\s+//g;
	s/>//;
	my ($chr, $element_start, $element_end, $element_length, $strand, $TE_ID, $TE_class);
	($chr, $element_start, $element_end, $TE_class) = ($1, $2, $3, $4) if /^(.*):([0-9]+)\.\.([0-9]+)#(.*)$/;
	$TE_ID = "$chr:$element_start..$element_end";
	$element_length = $element_end - $element_start + 1;
	next if $element_length < 80;
	$strand = ".";
	my $info = "ID=TE_annot_$i;Sequence_ontology=$SO{$class{$TE_class}};Name=$TE_ID#$TE_class;Method=structural";
        print GFF "$chr\t$annotator\t$TE_class\t$element_start\t$element_end\t.\t$strand\t.\t$info\n";
	$i++;
        }
close ID;
close GFF;

