#!/usr/bin/perl -w
use strict;

my $code = $ARGV[0]; #1, code name with serial numbers; 0, transform  back to original names with $seq.list
my $seq = $ARGV[1];

open Seq, "<$seq" or die $!;

$/ = "\n>";
my $i = 0;
my @seq;
my %seq;
while (<Seq>){
	s/>//g;
	chomp;
	my ($id, $seq) = (split /\n/, $_, 2);
	$seq =~ s/\s+//g;
	$seq[$i] = [$id, $seq];
	$seq{$id} = $seq;
	$i++;
	}
$/ = "\n";

if ($code == 1){
	open List, ">$seq.code.list";
	open Out, ">$seq.code";
	my $j = 0;
	foreach (@seq){
		print List "$j\t@{$_}[0]\n";
		print Out ">$j\n@{$_}[1]\n";
		$j++;
		}
	close List;
	close Out;
	}
elsif ($code == 0){
	open Code, "<$seq.list";
	open Out, ">$seq.decode";
	my %code;
	while (<Code>){
		chomp;
		my ($alt, $ori) = (split /\s+/, $_, 2);
		print Out ">$ori\n$seq{$alt}\n" if exists $seq{$alt};
		}
	close Code;
	close Out;
	}
