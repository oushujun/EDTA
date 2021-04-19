#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;
use File::Basename;
#Shujun Ou (shujun.ou.1@gmail.com)
#01/25/2021
#12/22/2019

my $usage = "\n
	Reclassify sequence based on RepeatMasker .out file
	The RM.out file is generated using a library with family classification to mask the seq.fa file.
		perl classify_by_lib_RM.pl -seq seq.fa -RM seq.fa.out -cov 80 -len 80 -iden 80\n\n";

#reclassify based on the 80-80-80 rule
my $min_iden = 80; #80%
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
open List, ">$seq.rename.list" or die $!;

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

	# check classifications
#	my ($sso, $qso) = ($class{"$s_class/$s_fam"}, "$q_class/$q_fam");
	my ($sso, $qso) = ("$s_class/$s_fam", $class{"$q_class/$q_fam"});
	if (exists $class{$sso}){
		$sso = $class{$sso}
		} else {
		print STDERR "$sso not found in the TE_SO database, it will not be used to rename sequences in the final annotation.\n";
		}
	next if $q_class ne $s_class and $q_class ne "NA";
	next if $qso ne $sso and $q_class ne "LTR";
#	next if "$q_class/$q_fam" ne "$s_class/$s_fam" and $q_class ne "NA";
	next if $q_class eq "DNA" and $q_fam =~ /Helitron|DHH/i and $s_fam !~ /Helitron|DHH/i;
	next if 100 - $div < $min_iden;
	my $len = $qe - $qs + 1;
	next if $len < $min_len;
	$subject =~ s/_INT|_LTR$//i; #combine LTR region and internal reigon into one family
	# note: if $qso=copia and $sso=other LTRs (gypsy or unknown), or $qso=gyp and $sso=others, or $qso=unk and $sso=others, they are all counted in coverage and thus the final $qso may not equal to $sso. The final fix will be improving the library and have these discripencies blocked.
	if (defined $lib{$query}){
		$lib{$query}{$subject} += $len;
		} else {
		$lib{$query} = {$subject => $len};
		}
	}
close RM;

foreach my $id (@lib){
	my $ori_id = $id;
	if (defined $lib{$id}){
		my $query_len = $lib{$id}{'len'};
		$lib{$id}{'len'} = 0;
		my @subjects = sort{$lib{$id}{$b} <=> $lib{$id}{$a}} (keys %{$lib{$id}});
		my ($topcov, $totcov) = (0, 0);
		$topcov = $lib{$id}{$subjects[0]}; #coverage of the query by the longest subject hit (%)
		$topcov = $topcov*100/$query_len;

		$totcov += $lib{$id}{$_} foreach @subjects;
		$totcov = $totcov*100/$query_len; #total coverage of the query by all subject hits (%)

		my $q_class = "NA";
		$q_class = $3 if $id =~ /:([0-9]+)\.\.([0-9]+)#[a-z]+\/([a-z]+)$/i;
		if ($q_class eq "Helitron" and @subjects > 1){ #rename this disregard coverage if it's a helitron
			$id = "$subjects[0]";
			}
		elsif ($totcov >= $min_cov and $topcov >= 30){ #inclusive parameter
#		elsif ($topcov >= $min_cov){
			$id = "$subjects[0]";
			}
		}
	print List "$ori_id\t$id\n";
	print New ">$id|$ori_id\n$seq{$ori_id}\n";
	}
close List;
close New;

