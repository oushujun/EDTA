##Shujun Ou
##usage: $ perl output_by_list.pl target_file list_file > outfile


#!usr/bin/perl -w
use strict;

my $usage="\n#usage: \$ perl output_by_list.pl DB_index_pos database LS_index_pos LIST [Exclusive]* [MSU_format] [FASTA_format] [version]> outfile
	* [] parameters are optional. 
		[Exclusive] -ex means exclude the entries in list, default is output the entries in list. 
		[MSU_format] -MSU0 means MSU_LOC occurs in the list file, while -MSU1 means MSU_LOC occurs in the database file
		eg. perl output_by_list.pl 1 Chr1.ltrTE.RMlist 1 Chr1.ltrTE.true.list -MSU0 -FA > Chr1.ltrTE.true.RMlist\n";
my $version="
output_by_list.pl
output_by_list: program for extracting information in database by provided list
Author: Shujun Ou, Department of Horticulture, Michigan State University, East Lansing, MI, 48823, USA
Version: 1.5 2014/05/12
\n";

my $msuL=0;
my $msuD=0; #1 for database file contain msu-type ID
my $exclude=0; #0 will output entries in the list from database; 1 will output entries not in the list from database
my $newline="\n"; #"\n>";

foreach my $para (@ARGV) {
	$msuL=1 if $para=~/^-MSU0$/i;
	$msuD=1 if $para=~/^-MSU1$/i;
	$exclude=1 if $para=~/^-ex$/i;
	$newline="\n>" if $para=~/^-FA$/i;
	die $version if $para=~/^-v$/i;
	die $usage if $para=~/^-u$/i;
	}

my $data_pos=$ARGV[0]-1;
open(DB,"$ARGV[1]") or die "ERROR: $!";
my $list_pos=$ARGV[2]-1;
open(LS,"$ARGV[3]") or die "ERROR: $!";


my %data;
while(<LS>){
	s/>//g;
	next if /^\s+$/;
	s/^\s+//;
	chomp;
	my $loc=(split)[$list_pos];
	$loc=~s/\|.*$//;
	$loc=~s/\[.*\]//g;
	if ($msuL){ #for MSU LOC position recognision 
		$loc=(split /:/, $loc)[1];
		$loc=~s/\.\..*$//;
		}
	$data{$loc}=undef;
}
close LS ;

$/="$newline";
while(<DB>){
	s/>//g;
	s/^\s+//;
	my $pos=(split)[$data_pos];
	if ($pos=~/pos/i){ print $_ }
	$pos=~s/\[.*\]//g;
	if ($pos=~/^([0-9]+),.*$/){$pos=$1}
	my ($p1, $p2)=($1, $2) if $pos=~/(.*)\|(.*)$/;
	if ($msuD){
		$pos=(split /:/, $pos)[1];
		$pos=~s/\.\..*$//;
		}
        if (exists $data{$pos} or exists $data{$p1} or exists $data{$p2}){
		if ($exclude==0){
			print ">" if $newline eq "\n>";
			print "$_";
			#sequences with the same ID will be output only once
			delete $data{$pos};
			delete $data{$p1};
			delete $data{$p2};
			}
		} else {
		if ($exclude==1){
			print ">" if $newline eq "\n>";
			print "$_";
			}
		}
}
close DB;


