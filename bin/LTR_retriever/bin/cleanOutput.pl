#!/usr/bin/env perl -w
use strict;

my ($index, $LINE, $DNA, $PlantP)=('', '', '', '');
($index, $LINE, $DNA, $PlantP)=@ARGV[0,1,2,3];
die "ERROR: Usage: perl cleanOut.pl index\n" if $index eq '';

`rm $LINE* $DNA* $PlantP* 2>/dev/null`;

`rm $index.cat.gz $index.LTRlib.clust.clstr $index.LTRlib $index.LTRlib.raw $index.LTRID.list $index.LTRlib.fa.n* $index.ltrTE* $index.nmtf $index.prelib* $index.nmtf.prelib $index.retriever.scn* $index.retriever.all.scn.list $index.*exclude* $index.*.nhr $index.*.nin $index.*.nsq 2>/dev/null`;

