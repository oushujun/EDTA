#!/usr/bin/env perl
use strict;
use warnings;

# usage: perl count_base.pl TElib.fa -s | this_script.pl > TElib.fa.info
# Shujun Ou (07/04/2021) shujun.ou.1@gmail.com

my %info;
while (<>){
	next if /^All\s+/;
	my ($id, $len, $miss) = (split)[0,1,2];
	$id =~ s/(_INT|_LTR)?#.*//;
	$len -= $miss;
	$info{$id} += $len;
	}

print "id\tfam_len\n";
foreach my $id (keys %info){
	print "$id\t$info{$id}\n";
	}
