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
my $num = 0;
$num = $ARGV[1] if defined $ARGV[1] and $ARGV[1] =~ /^[0-9]+$/;
my %data;
my @map; # store TE_name => original_name mappings
while (<FA>){
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
		if (exists $data{$loc}){
			my $record_num = (split /\s+/, $data{$loc})[0];
			$record_num =~ s/>TE_([0-9]+)_(LTR|INT).*/$1/i;
			$data{$loc} .= ">TE_${record_num}_$part#$class\n$seq\n";
			#print "add: $loc\t$record_num\t$part#$class\n"; #test
			} else {
			$data{$loc} = ">TE_${num}_$part#$class\n$seq\n";
			push @map, "TE_${num}_$part#$class\t$name";
			#print "new: $loc\t$num\t$part#$class\n"; #test
			$num++;
			}
		} else {
		$data{$num} = ">TE_${num}#$class\n$seq\n";
		push @map, "TE_${num}#$class\t$name";
		$num++;
		}
	}
close FA;

foreach my $fam (sort{$data{$a} cmp $data{$b}} (keys %data)){
	print $data{$fam};
	}

# write mapping file if requested
if ($mapfile ne ''){
	open MAP, ">$mapfile" or die "Cannot open $mapfile for writing: $!\n";
	print MAP "$_\n" for @map;
	close MAP;
}