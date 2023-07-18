#!/usr/bin/env perl
use warnings;
use strict;
use File::Basename;

my $usage = "\n\tperl bed2gff.pl EDTA.TE.combo.bed > EDTA.TE.combo.gff3\n\n";
my $annotator="EDTA";

#obtain the exact path for the program location
my $script_path = dirname(__FILE__);

my $ID_name = "TE_annot";
$ID_name = $ARGV[1] if defined $ARGV[1]; # TE_homo or TE_struc

#read SO info and defined sequence ontology
my $SO = "$script_path/TE_Sequence_Ontology.txt";
open SO, "<$SO" or die "The sequence ontology file 'TE_Sequence_Ontology.txt' is not found in $script_path!\n";
my (%class, %SO);
while (<SO>){
	next if /#/;
	next if /^(\s+)?$/;
	my ($so_name, $so_id, $so_alias) = (split /\s+/, $_, 3);
	$so_alias =~ s/\s+//;
	$SO{$so_name} = $so_id;
	$class{$so_name} = $so_name;
	foreach my $alia ((split /,/, $so_alias)){
		$class{$alia} = $so_name;
		}
	}
close SO;

my $date=`date -u`;
chomp ($date);
open BED, "sort -sV -k1,1 -k2,2 $ARGV[0] |" or die $usage;
print "##gff-version 3\n##date $date
##seqid source sequence_ontology start end score strand phase attributes\n";

my $i = 0; #annotation order
while (<BED>){
	next if /^#/;
	chomp;
	my ($chr, $element_start, $element_end, $element_length, $TE_class, $method, $iden, $score, $strand, $phase, $TE_ID, $info, $extra);
	($chr, $element_start, $element_end, $TE_ID, $TE_class, $method, $iden, $score, $strand, $phase, $extra) = (split);
	$element_length = $element_end - $element_start + 1;

	my $so = $TE_class;
	if (exists $class{$TE_class}){
		$so = $class{$TE_class}
		} else {
		print STDERR "\nWarning: $TE_class not found in the TE_SO database, will use the general term 'repeat_region\tSO:0000657' to replace it.\n";
		$so = "repeat_region";
		}
	$extra = '' unless defined $extra;
	$extra = '' if $extra eq 'NA' or $extra eq '.';
	if ($extra eq ''){
		$extra = '';
		} else {
		$extra = ";$extra";
		}
	$info = "ID=${ID_name}_$i;Name=$TE_ID;Classification=$TE_class;Sequence_ontology=$SO{$so};Identity=$iden;Method=${method}$extra";
	print "$chr\t$annotator\t$so\t$element_start\t$element_end\t$score\t$strand\t$phase\t$info\n";
	$i++;
	}
close BED;
