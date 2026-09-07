#!/usr/bin/env perl
use warnings;
use strict;

my $usage = "Usage: perl rename_TE.pl input.fa [start_num] [--map mapfile]\n";
my $fasta = $ARGV[0];
my $mapfile = '';

# parse arguments
my $k = 0;
foreach (@ARGV){
	$mapfile = $ARGV[$k+1] if /^--map$/i and defined $ARGV[$k+1];
	$k++;
}

open FA, "<$fasta" or die "\nInput not found!\n$usage";
$/ = "\n>";
my @records = <FA>;
close FA;
# canonical input order: upstream pools (e.g. the RepeatModeler-derived novel
# sequences) do not guarantee a stable record order, and TE numbering depends on
# input order — sort by header so identical sequence sets always number identically.
# strip ">" from the keys: only the first record keeps one after the "\n>" split
@records = sort { my ($x) = split /\n/, $a; my ($y) = split /\n/, $b; $x =~ s/>//g; $y =~ s/>//g; $x cmp $y or $a cmp $b } @records;
my $num = 0;
$num = $ARGV[1] if defined $ARGV[1] and $ARGV[1] =~ /^[0-9]+$/;
my %data;
my %pair; #$loc => parts stored under the current element number, so reappearing parts get a new number
my %curr_num; #$loc => the element number currently in use
my %map_lines; #output group => map lines, filled in output order
foreach my $record (@records){
	$_ = $record;
	s/>//g;
	$num = sprintf("%08d", $num);
	my ($id, $seq) = (split /\n/, $_, 2);
	my $name = (split /\s+/, $id)[0];
	$seq =~ s/\s+//g;
	my ($fam, $class) = ($1, $2) if $name =~ /^(.*)#(.*)$/;
	#print "$id\t$fam, $class\n"; #test
	#rename TE as unknown if $class info could not be retrieved
	$class = "unknown" unless defined $class;
	#retain LTR-INT info for LTR sequences
	if ($name =~ /_(LTR|INT)#LTR/i){
		my ($loc, $part) = ('', '');
		($loc, $part) = ($1, $2) if $fam =~ /^(.*)_(LTR|INT)$/i;
		$loc = "$loc#$class";
		#print "Ori: $name\t$loc\n"; #test
		if (not exists $curr_num{$loc} or exists $pair{$loc}{$part}){
			#this family+part has reappeared: allocate a unique number for the new occurrence
			$curr_num{$loc} = $num;
			$pair{$loc} = {};
			$num++;
			}
		$pair{$loc}{$part} = 1;
		$data{$loc} .= ">TE_$curr_num{$loc}_$part#$class\n$seq\n";
		push @{ $map_lines{$loc} }, "TE_$curr_num{$loc}_$part#$class\t$name";
		#print "add: $loc\t$curr_num{$loc}\t$part#$class\n"; #test
		} else {
		$data{$num} = ">TE_${num}#$class\n$seq\n";
		push @{ $map_lines{$num} }, "TE_${num}#$class\t$name";
		$num++;
		}
	}

my @map;
foreach my $fam (sort{$data{$a} cmp $data{$b}} (keys %data)){
	print $data{$fam};
	push @map, @{ $map_lines{$fam} };
	}

# write mapping file if requested
if ($mapfile ne ''){
	open MAP, ">$mapfile" or die "Cannot open $mapfile for writing: $!\n";
	print MAP "$_\n" for @map;
	close MAP;
}