#!/usr/bin/env perl
use warnings;
use strict;

# Convert gff files to enriched bed format
# Shujun Ou (shujun.ou.1@gmail.com)
# v0.6	12/21/2023
# v0.5	04/21/2023
# v0.4	09/27/2021
# v0.3	07/22/2020
# v0.2	07/20/2020
# v0.1	12/19/2019

my $usage = "Usage:
	perl gff2bed.pl file.gff [structural|homology] > file.bed
	cat *.gff3 | perl gff2bed.pl - [structural|homology] > file.bed \n";
die "Please indicate if the specified GFF file is generated based on structural features or homology!\n$usage\n" unless defined $ARGV[1] and $ARGV[1] =~ /^structural$|^homology$/i;
my $method = $ARGV[1];
open GFF, "<$ARGV[0]" or die $usage;

while (<GFF>){
	chmod;
	next if /^#/;

	my ($method, $type, $TE_class, $class, $iden) = (undef, undef, undef, undef, 'NA');
	my ($chr, $sequence_ontology, $element_start, $element_end, $score, $strand, $phase, $extra) = (split)[0,2,3,4,5,6,7,8];
	next unless defined $chr and defined $element_start;

	# get class info for summary categories
	$class = $sequence_ontology;
	$class = 'non_LTRretrotransposon' if $class eq "non_LTR_retrotransposon"; #avoid confusing the LTR processing
	# combine LTR subclass
	if ($class =~ /(LTR_retrotransposon|_LTR_retrotransposon|long_terminal_repeat|TRIM|LARD)/i){
		$class =~ s/_LTR_retrotransposon//i;
		$class = "unknown" if $class eq "LTR_retrotransposon";
		$class = "LTR/$class";
		$type = "LTR";
		}
	# combine TIR subclass
	if ($class =~ /(terminal_inverted_repeat|TIR_transposon|polinton|MITE)/i){
		$class =~ s/_TIR_transposon//i;
		$class =~ s/terminal_inverted_repeat_element/unknown/i;
		$class =~ s/terminal_inverted_repeat/TIR/i;
		$class = "TIR/$class";
		$type = "TIR";
		}
	# combine LINE subclass
	if ($class =~ /(LINE|LINE_element|LINE_retrotransposon|RIL)/i){
		$class = "unknown" if $class eq "LINE_element";
		$class =~ s/_LINE_retrotransposon//i;
		$class = "LINE/$class";
		$type = "LINE";
		}
	# combine SINE subclass
	if ($class =~ /(SINE|SINE_element|SINE_retrotransposon|RIS)/i){
		$class = "unknown" if $class eq "SINE_element";
		$class =~ s/_SINE_retrotransposon//i;
		$class = "SINE/$class";
		$type = "SINE";
		}
	if ($class =~ /(nonLTR|non_LTR|YR_retrotransposon|Penelope|ERTBV|pararetrovirus)/i){
		$class = "unknown" if $class eq "non_LTRretrotransposon";
		$class =~ s/_retrotransposon//i;
		$class = "nonLTR/$class";
		$type = "nonLTR";
		}
	if ($class =~ /(rDNA|rRNA|rDNA_intergenic_spacer_element|rRNA_gene)/i){
		$class = "rDNA/IGS" if $class eq "rDNA_intergenic_spacer_element";
		$class = "rDNA/45S" if $class eq "rRNA_gene";
		$class =~ s/cytosolic_(.*)_rRNA/rDNA\/$1/i;
		$class =~ s/rRNA_([0-9])_external_transcribed_sequence/rDNA\/$1_ETS/i;
		$class =~ s/rRNA_internal_transcribed_spacer([0-9])/rDNA\/ITS$1/i;
		$type = "rDNA";
		}
	$class = "nonTIR/$class" if $class =~ /(Crypton_YR_transposon|helitron)/i;

	# determine $type for struc-homo TE annotation merging
	# $type is critical for get_frag.pl and keep_nest.pl, 
	# if this info is not properlly assigned, the resulting bed file is wrong and EDTA.TEanno.sum will be empty
	# $type should match what's listed in util/TE_Sequence_Ontology.txt
	$type = "Cent" if $sequence_ontology =~ /Cent|centromeric_repeat/i;
	$type = "knob" if $sequence_ontology =~ /knob/i;
	$type = "satellite" if $sequence_ontology =~ /satellite|satellite_DNA/i;
	$type = "telomere" if $sequence_ontology =~ /telomer|telomeric_repeat/i;
	$type = "subtelomere" if $sequence_ontology =~ /subtelomer/i;
	$type = "low_complexity" if $sequence_ontology =~ /low_complexity/i;
	$type = "Helitron" if $sequence_ontology =~ /Helitron|DHH/i;
	$type = "Crypton" if $sequence_ontology =~ /Crypton_YR_transposon/i;
	$type = "plastid" if $sequence_ontology =~ /chloroplast|plastid|mitochondri/i;
	$type = "repeat_region" if $sequence_ontology =~ /(repeat_region|DNA_transposon|^retrotransposon$)/i;
	$type = 'repeat_region' if $sequence_ontology =~ /Unknown|Unspecified/i; #suggested by Changfu Jia
	$type = $1 if $sequence_ontology =~ /^(.*)\/.*/ and $1 !~ /DNA|MITE/i;

	# get assortive structural info
	my $TE_ID = "$chr:$element_start..$element_end";
	$TE_ID = $1 if $extra =~ s/Name=(.*?);//i;
	$TE_class = $1 if $extra =~ s/Classification=(.*?);//i;
	$iden = $1 if $extra =~ s/ltr_identity=([0-9.e\-]+);//i or $extra =~ s/Identity=([0-9.e\-]+);//i;
	$method = $1 if $extra =~ s/Method=(homology|structural)//i;
	$extra =~ s/ID=.*Sequence_ontology=SO:[0-9]+;//; #rename annotation id based on input order
	$extra =~ s/^;//;
	$extra =~ s/;$//;
	$extra =~ s/;+/;/g;
	$extra = "NA" if $extra =~/^$/;

	# skip some structural TE entries
	next if $sequence_ontology =~ /(target_site_duplication|primer_binding_site|U_box|RR_tract)/i;
	next if $sequence_ontology eq "repeat_region" and $TE_class =~ /LTR/i;
	next if $sequence_ontology eq "long_terminal_repeat" and $method =~ /structural/i;
	next if $sequence_ontology eq "terminal_inverted_repeat" and $method =~ /structural/i;

	print "$chr\t$element_start\t$element_end\t$TE_ID\t$TE_class\t$method\t$iden\t$score\t$strand\t$phase\t$extra\t$type\t$class\n";
	}

