#!/usr/bin/env perl
use warnings;
use strict;

my $usage = "";
my $fasta = $ARGV[0];

open FA, "<$fasta" or die "\nInput not found!\n$usage";
$/ = "\n>";
my $num = 0;
$num = $ARGV[1] if defined $ARGV[1];
my %data;
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
			#print "new: $loc\t$num\t$part#$class\n"; #test
			$num++;
			}
		} else {
		$data{$num} = ">TE_${num}#$class\n$seq\n";
		$num++;
		}
	}
close FA;

foreach my $fam (sort{$data{$a} cmp $data{$b}} (keys %data)){
	print $data{$fam};
	}

