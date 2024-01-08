#!/usr/bin/env perl
use warnings;
use strict;
# Shujun Ou (shujun.ou.1@gmail.com)
# 12/24/2019 @ DSM airport

my $usage = "
	Replace strings in the target with a conversion list (old, new)
		perl rename_by_list.pl target list mode[0|1]
			mode = 0, generic replace (slow)
			mode = 1, specific for gff3 files\n\n";

my $target = $ARGV[0];
my $list = $ARGV[1];
my $mode = 0;
$mode = $ARGV[2] if defined $ARGV[2];

open List, "<$list" or die $usage;
open Target, "<$target" or die $usage;

my %list;
while (<List>){
	chomp;
	next if /^#/;
	my ($old_id, $new_id) = (split)[0,1];
	$list{$old_id} = $new_id;
	$old_id =~ s/#.*//; #make old_id flexible by removing the classification info
	$new_id =~ s/#.*//;
	$list{$old_id} = $new_id;
	}
close List;

while (<Target>){
	if ($mode eq 0){
		foreach my $old_id (keys %list){
			my $new_id = $list{$old_id};
			s/$old_id/$new_id/g;
			}
		}
	elsif ($mode eq 1){
		my $id = $2 if /(Name)=(.*?);/i;
		my $new = $list{$id} if defined $id and defined $list{$id};
		s/$id/$new/g if defined $new;
		}
	print $_;
	}
close Target;
