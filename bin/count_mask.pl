#!/usr/bin/env perl
use warnings;
use strict;

my $add = 0;
while (<>){
	my ($s,$e)=(split)[1,2];
	my $len=$e-$s+1;
	$add+=$len;
	}
print "$add\n";
