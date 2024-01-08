#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;
use File::Basename;

#Shujun Ou (shujun.ou.1@gmail.com)
#01/06/2024 facilitated by ChatGPT
#01/25/2021
#12/22/2019

my $usage = "\n
	Reclassify sequence based on RepeatMasker .out file
	If new classification is not the same as the original, the family name won't be changed and will be exported into the false list.
	The RM.out file is generated using a library with family classification to mask the seq.fa file.
		perl classify_by_lib_RM.pl -seq seq.fa -RM seq.fa.out -cov 80 -len 80 -iden 80\n\n";

#reclassify based on the 80-80-60 rule
my $min_iden = 60; #60%, sufficient for RepeatMasker alignments, which usually are the best hits already
my $min_len = 80; #80bp
my $min_cov = 80; #80%

my $seq = '';
my $RM = '';

my $k = 0;
foreach (@ARGV){
	$seq = $ARGV[$k+1] if /^-seq$/i;
	$RM = $ARGV[$k+1] if /^-RM$/i;
	$min_cov = $ARGV[$k+1] if /^-cov$/i;
	$min_len = $ARGV[$k+1] if /^-len$/i;
	$min_iden = $ARGV[$k+1] if /^-iden$/i;
	die $usage if /^-h$|^-help$|^--help$/i;
	$k++;
	}

open Seq, "<$seq" or die $usage;
open RM, "<$RM" or die $usage;
open New, ">$seq.rename" or die $!;
open False, ">$seq.false.list" or die $!;
open Good, ">$seq.rename.list" or die $!;
open List, ">$seq.RMall.list" or die $!;
print List "#Original\tNew\tQuery_ontology\tSubject_ontology\tTopHit_cov\tTotal_cov\n";
print Good "#Original\tNew\tQuery_ontology\tSubject_ontology\tTopHit_cov\tTotal_cov\n";
print False "#Original\tNew\tQuery_ontology\tSubject_ontology\tTopHit_cov\tTotal_cov\n";

#obtain the exact path for the program location
my $script_path = dirname(__FILE__);

##read SO info and defined sequence ontology
my $SO = "$script_path/TE_Sequence_Ontology.txt";
open SO, "<$SO" or die "The sequence ontology file 'TE_Sequence_Ontology.txt' is not found in $script_path!\n";
my (%class, %SO);
while (<SO>){
	next if /#/;
	next if /^(\s+)?$/;
	my ($so_name, $so_id, $so_alias) = (split /\s+/, $_, 3);
	$so_alias =~ s/\s+//;
	$SO{$so_name} = $so_id;
	foreach my $alia ((split /,/, $so_alias)){
		$class{$alia} = $so_name;
		}
	}
close SO;

# read in lib seqs
my %lib;
my %seq;
my @lib;
$/ = "\n>";
while (<Seq>){
	s/>//g;
	my ($id, $seq) = (split /\n/, $_, 2);
	$id = (split /\s+/, $id)[0];
	$seq =~ s/\s+//g;
	$seq{$id} = $seq;
	my $len = length $seq;
	$lib{$id} = {'len' => $len};
	push @lib, $id;
	}
close Seq;
$/ = "\n";

# process RM.out data
while (<RM>){
	next if /^\s+?$/;
	my @hit = (split);
	my ($query, $subject, $div, $qs, $qe, $s_class) = @hit[4,9,1,5,6,10];
	next unless exists $lib{$query};
	my ($q_class, $q_fam, $s_fam) = ("NA", "NA", "NA");
	($s_class, $s_fam) = ($1, $2) if $s_class =~ /^(.*)\/(.*)$/;
	($q_class, $q_fam) = ($1, $2) if $query =~ /#(.*)\/(.*)$/;

	# standardize query and subject classifications using the TE_SO database
	my ($sso, $qso) = ("$s_class/$s_fam", $class{"$q_class/$q_fam"});
	if (exists $class{$sso}){
		$sso = $class{$sso}
		} else {
		print STDERR "$sso not found in the TE_SO database, it will not be used to rename sequences in the final annotation.\n";
		}
	next if 100 - $div < $min_iden;
	my $len = $qe - $qs + 1;
	next if $len < $min_len;

	# combine LTR region and internal reigon into one family
	# note: if $qso=copia and $sso=other LTRs (gypsy or unknown), or $qso=gyp and $sso=others, or $qso=unk and $sso=others, they are all counted in coverage and thus the final $qso may not equal to $sso. The final fix will be improving the library and have these discripencies blocked.
	$subject =~ s/(_I_nonauto|_I_auto|-I_nonauto|-I_nona|-I_auto|_LTR_nonauto|_LTR_auto|\-LTR|\-I|_I|_INT|_LTR|I|LTR)$//i;

	# store query-sso-subject matches. SSO is the bigger category including many subjects. For each query-subject match, here stores the matched length. 
	if (defined $lib{$query}{$sso}{$subject}){
		$lib{$query}{$sso}{$subject} += $len;
		} else {
		$lib{$query}{$sso}{$subject} = $len;
		}
	}
close RM;

# Find out the sso group that accumulates the longest match to query, and which subject under this sso provides the longest match. 
foreach my $id (@lib){
	# obtain sequence ontology of this query id (qso)
	my ($ori_id, $q_class, $qso) = ($id, "NA", "NA");
	$q_class = $1 if $id =~ /#(.*)$/i;
	$qso = $class{"$q_class"} if defined $class{"$q_class"};
	my ($query_len, $top_coverage, $total_coverage) = ($lib{$id}{'len'}, 0, 0);
	my ($total_subject_length, $max_sso, $max_sso_length, $max_subject, $max_subject_length) = (0, "NA", 0, "NA", 0);
	if (defined $lib{$id}){
		# store sso and subject with their lengths
		my %sso_subjects = %{ $lib{$id} };
		foreach my $sso (keys %sso_subjects) {
			next if $sso eq 'len';
			my $total_sso_length = 0;
			my $max_length_subject = "";
			my $max_length = 0;
			foreach my $subject (keys %{ $sso_subjects{$sso} }) {
				my $length = $sso_subjects{$sso}{$subject};
				$total_sso_length += $length;
				$total_subject_length += $length;
				if ($length > $max_length) {
					$max_length = $length;
					$max_length_subject = $subject;
					}
				}

			if ($total_sso_length > $max_sso_length) {
				$max_sso_length = $total_sso_length;
				$max_sso = $sso;
				$max_subject = $max_length_subject;
				$max_subject_length = $max_length;
				}
			}

		#my @subjects = sort{$lib{$id}{$b} <=> $lib{$id}{$a}} (keys %{$lib{$id}});
		$top_coverage = sprintf("%.1f", $max_subject_length*100/$query_len); #coverage of the query by the longest subject hit (%)
		$total_coverage = sprintf("%.1f", $total_subject_length*100/$query_len); #total coverage of the query by all subject hits (%)
		
		# rename $id to the longest alignment if passing criteria
		if ($q_class =~ /Helitron/i and keys %sso_subjects > 1){ #rename this disregard coverage if it's a helitron
			$id = $max_subject;
			}
		}
	
	# remove classification in the new id, store all matched and nonmatched cases
	$id =~ s/#.*//;
	print List "$ori_id\t$id\t$qso\t$max_sso\t$top_coverage\t$total_coverage\n";

	# decide when to take the new classification, use inclusive parameters
	if ($max_sso ne "NA" and $total_coverage >= $min_cov and $top_coverage >= 30){
		# new classification match the old classification, thus rename
		$id = $max_subject;
		if (($qso eq $max_sso) or 
		    ($qso =~ /LTR_retrotransposon/ && $max_sso =~ /LTR_retrotransposon/) or
		    ($qso =~ /TIR_transposon/ && $max_sso =~ /TIR_transposon/) or
		    ($qso =~ /LINE_/ && $max_sso =~ /LINE_/) or
		    ($qso =~ /SINE_/ && $max_sso =~ /SINE_/) or
		    ($qso =~ /rRNA/ && $max_sso =~ /rRNA/) or
		    ($qso =~ /rDNA/ && $max_sso =~ /rDNA/) or
		    ($qso =~ /rRNA/ && $max_sso =~ /rDNA/) or
		    ($qso =~ /rDNA/ && $max_sso =~ /rRNA/) or
		    ($qso =~ /YR_retrotransposon/ && $max_sso =~ /YR_retrotransposon/)){
			print Good "$ori_id\t$id\t$qso\t$max_sso\t$top_coverage\t$total_coverage\n";
			} else {
			# new classification does NOT match the old classification, thus reverse to ori_id
			$id = $1 if $ori_id =~ /^(.*)#/;
			print False "$ori_id\t$id\t$qso\t$max_sso\t$top_coverage\t$total_coverage\n";
			}

		} else {
		# no match or below threshold, thus reverse to ori_id
		$id = $1 if $ori_id =~ /^(.*)#/;
		}
	#print New ">$id#$q_class\t$ori_id\n$seq{$ori_id}\n";
	print New ">$id|$ori_id\n$seq{$ori_id}\n";
	}
close List;
close Good;
close False;
close New;

