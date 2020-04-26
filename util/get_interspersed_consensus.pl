#!/usr/bin/env perl
use warnings;
use strict;

my $file = "interspersed_repeat.alignment";
my $minlen = 80;

open File, "<$file" or die $!;

my $new = 0;
my $count = 0;
while (<File>){
	chomp;
	($new = 1 and next) if /^\-+/;
	($new = 0 and next ) if /^>/;
	next if length $_ < $minlen;
	print ">$count\n$_\n" and $count++ if $new == 1;
	}
close File;
