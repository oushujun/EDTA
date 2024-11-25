#!/usr/bin/env perl
use warnings;
use strict;
use File::Basename;

# Format sequence IDs generated from other programs using the TE SO system.
# Shujun Ou (shujun.ou.1@gmail.com)
# v0.1	06/24/2021

my $usage = "\n\tperl format_TElib.pl EDTA.TElib.fa > EDTA.TElib.fa.mod\n\n";

my $allow_unknown = 1; #when the TE_SO does not have info: 1, keep the original id; 0, replace to "repeat_region".

#obtain the exact path for the program location
my $script_path = dirname(__FILE__);

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

# process the fasta file
open FA, "<$ARGV[0]" or die $usage;
$/="\n>";
while (<FA>){
	chomp;
	s/>//g;
	my ($id, $seq)=(split /\n/, $_, 2);
	$seq=~s/\s+//g;
	my $TE_class = $1 if $id =~ s/#(.*)$//;
	my $supfam = $TE_class;
	$supfam = 'repeat_region/Unknown' if $supfam =~ /^unknown$/i;
	if ($TE_class =~ /\/unknown$/i){
		$supfam = $1 if $TE_class =~ /^(.*)\//;
		} else {
		$supfam = $1 if $TE_class =~ s/^(.*)\///;
		}
	$supfam = 'TIR' if $supfam =~ /^DNA$|^MITE$/i and $TE_class !~ /Helitron/i;

	my $so = $TE_class;
	if (exists $class{$TE_class}){
		$so = $class{$TE_class}
		} 
	elsif ($allow_unknown == 0) {
		print STDERR "\nWarning: $TE_class not found in the TE_SO database, will use the general term 'repeat_region\tSO:0000657' to replace it.\n";
		$so = "repeat_region";
		}
	$id = "$id#$supfam/$so";
	$id =~ s/Unknown\/repeat_region/Unknown/i;

	print ">$id\n$seq\n";
	}
close FA;
