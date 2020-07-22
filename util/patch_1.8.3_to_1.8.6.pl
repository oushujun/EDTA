#!/usr/bin/env perl
use strict;
use warnings;
use FindBin;

my $genome = $ARGV[0];
$genome = "$genome.mod";
my $threads = 4;
$threads = $ARGV[1] if defined $ARGV[1] and $ARGV[1] =~ /^\d+$/;

#dependencies
my $script_path = $FindBin::Bin;
my $filter_gff = "$script_path/filter_gff3.pl";
my $make_gff3 = "~/las/git_bin/LTR_retriever/bin/make_gff3.pl";
my $output_by_list = "$script_path/output_by_list.pl";
my $bed2gff = "$script_path/bed2gff.pl";
my $make_bed = "$script_path/make_bed_with_intact.pl";
my $RMout2bed = "$script_path/RMout2bed.pl";
my $rename_by_list = "$script_path/rename_by_list.pl";
my $gff2bed = "$script_path/gff2bed.pl";
my $get_frag = "$script_path/get_frag.pl";
my $keep_nest = "$script_path/keep_nest.pl";
my $combine_overlap = "$script_path/combine_overlap.pl";
my $count_base = "$script_path/count_base.pl";
my $buildSummary = "$script_path/util/buildSummary.pl";


if (0){
## LTR
chdir "$genome.EDTA.raw/LTR";

`perl $make_gff3 $genome $genome.pass.list`;
`perl $output_by_list 1 $genome.LTR.intact.fa.ori 1 $genome.LTR.intact.fa -FA -MSU0 -MSU1 -ex|grep \\>|perl -nle 's/>//; s/#.*//; print "Name\\t\$_"' > $genome.LTR.intact.fa.ori.rmlist`; #update
`perl $filter_gff $genome.pass.list.gff3 $genome.LTR.intact.fa.ori.rmlist | perl -nle 's/LTR_retriever/EDTA/gi; print \$_' > $genome.LTR.intact.gff3`;
`cp $genome.LTR.intact.gff3 ../`;
chdir "../../";

## TIR
chdir "$genome.EDTA.raw/TIR";

`perl -nle 's/\\-\\+\\-/_Len:/; my (\$chr, \$method, \$supfam, \$s, \$e, \$anno) = (split)[0,1,2,3,4,8]; my \$class='DNA'; \$class='MITE' if \$e-\$s+1 <= 600; my (\$tir, \$iden, \$tsd)=(\$1, \$2/100, \$3) if \$anno=~/TIR:(.*)_([0-9.]+)_TSD:([a-z0-9._]+)_LEN/i; print "\$chr \$s \$e \$chr:\$s..\$e \$class/\$supfam structural \$iden . . . TSD=\$tsd;TIR=\$tir"' ./TIR-Learner-Result/TIR-Learner_FinalAnn.gff3 | perl $output_by_list 4 - 1 $genome.TIR.raw.fa -MSU0 -MSU1 > $genome.TIR.intact.bed`;
`perl $bed2gff $genome.TIR.intact.bed TIR > $genome.TIR.intact.gff3`;
`cp $genome.TIR.intact.bed $genome.TIR.intact.gff3 ../`;
chdir "../../";

## Helitron
chdir "$genome.EDTA.raw/Helitron";

`perl -nle 's/^(>.*)\\s+(.*)\$/\$1#DNA\\/Helitron\\t\$2/; print \$_' $genome.HelitronScanner.filtered.fa | perl $output_by_list 1 - 1 $genome.Helitron.intact.fa -FA -MSU0 -MSU1 > $genome.Helitron.intact.fa.temp`;
`perl $make_bed $genome.Helitron.intact.fa.temp > $genome.Helitron.intact.bed`;
`perl $bed2gff $genome.Helitron.intact.bed HEL > $genome.Helitron.intact.gff3`;
`cp $genome.Helitron.intact.bed $genome.Helitron.intact.gff3 ../`;
chdir "../../";

## raw
chdir "$genome.EDTA.raw";

`cat $genome.TIR.intact.bed $genome.Helitron.intact.bed | perl $bed2gff - TE_struc > $genome.EDTA.intact.gff3.raw`;
`cat $genome.LTR.intact.gff3 >> $genome.EDTA.intact.gff3.raw`;
`sort -sV -k1,1 -k4,4 $genome.EDTA.intact.gff3.raw | grep -v '^#' > $genome.EDTA.intact.gff3; rm $genome.EDTA.intact.gff3.raw`;
`cp $genome.EDTA.intact.gff3 ../`;
chdir "../";

## final
chdir "$genome.EDTA.final";

`cp ../$genome.EDTA.raw/$genome.EDTA.intact.fa ./$genome.EDTA.intact.fa.raw`;
`cp ../$genome.EDTA.raw/$genome.EDTA.intact.gff3 ./`;

my $ls = `ls $genome*masked.cleanup.rmlist`;
if ($ls =~ /masked.cleanup.rmlist/){
	`cat $genome*masked.cleanup.rmlist | perl -nle 's/#.*//; print \$_' > $genome.EDTA.intact.fa.masked.cleanup.rmlist.all`; #update
	`perl $filter_gff $genome.EDTA.intact.gff3 $genome.EDTA.intact.fa.masked.cleanup.rmlist.all > $genome.EDTA.intact.gff3.new`;
	`perl -nle 'my \$id = \$1 if /=(repeat_region[0-9]+);/; print "Parent\t\$id\nName\t\$id" if defined \$id' $genome.EDTA.intact.gff3.removed >> $genome.EDTA.intact.fa.masked.cleanup.rmlist.all`;
	`perl $filter_gff $genome.EDTA.intact.gff3 $genome.EDTA.intact.fa.masked.cleanup.rmlist.all > $genome.EDTA.intact.gff3.new`;
	`mv $genome.EDTA.intact.gff3.new $genome.EDTA.intact.gff3`; #update intact.gff
	}

if (-s "$genome.EDTA.intact.fa.rename.list"){
        `perl $rename_by_list $genome.EDTA.intact.gff3 $genome.EDTA.intact.fa.rename.list 1 > $genome.EDTA.intact.gff3.rename`;
        `mv $genome.EDTA.intact.gff3.rename $genome.EDTA.intact.gff3`;
        `cp $genome.EDTA.intact.gff3 ../`; #replace the intact gff that has no lib family info
 	}
chdir "../";
}
#exit;


## Anno
chdir "$genome.EDTA.anno";
`ln -s ../$genome $genome` unless -e $genome;

if (0){
`cp ../$genome.EDTA.final/$genome.EDTA.intact.gff3 ./`;
`perl $RMout2bed $genome.EDTA.RM.out > $genome.EDTA.RM.bed`;
`perl $bed2gff $genome.EDTA.RM.bed TE_homo > $genome.EDTA.RM.gff3`;
`perl $gff2bed $genome.EDTA.intact.gff3 structural > $genome.EDTA.intact.bed`;
`perl $gff2bed $genome.EDTA.RM.gff3 homology > $genome.EDTA.RM.bed`;
`perl $combine_overlap $genome.EDTA.intact.bed $genome.EDTA.intact.bed.cmb 5`;
#`perl $combine_overlap $genome.EDTA.RM.bed $genome.EDTA.RM.bed.cmb 5`;
`perl $get_frag $genome.EDTA.RM.bed $genome.EDTA.intact.bed.cmb $threads`;
`perl $keep_nest $genome.EDTA.intact.bed $genome.EDTA.RM.bed $threads`;
`sort -suV $genome.EDTA.intact.bed-$genome.EDTA.RM.bed $genome.EDTA.RM.bed-$genome.EDTA.intact.bed.cmb > $genome.EDTA.TEanno.bed`;
`grep homology $genome.EDTA.intact.bed-$genome.EDTA.RM.bed > $genome.EDTA.intact.bed-$genome.EDTA.RM.bed.homo`;
`sort -suV $genome.EDTA.intact.bed-$genome.EDTA.RM.bed.homo $genome.EDTA.RM.bed-$genome.EDTA.intact.bed.cmb > $genome.EDTA.homo.bed`;
`perl $bed2gff $genome.EDTA.homo.bed TE_homo > $genome.EDTA.homo.gff3`;
`cat $genome.EDTA.intact.gff3 $genome.EDTA.homo.gff3 > $genome.EDTA.TEanno.gff3.raw`;
`grep -v '^#' $genome.EDTA.TEanno.gff3.raw | sort -sV -k1,1 -k4,4 | perl -0777 -ne '\$date=\`date\`; \$date=~s/\\s+\$//; print "##gff-version 3\\n##date \$date\\n##Identity: Sequence identity (0-1) between the library sequence and the target region.\\n##ltr_identity: Sequence identity (0-1) between the left and right LTR regions.\\n##tsd: target site duplication.\\n##seqid source sequence_ontology start end score strand phase attributes\\n\$_"' - > $genome.EDTA.TEanno.gff3`;
`rm $genome.EDTA.TEanno.gff3.raw`;
`cp $genome.EDTA.TEanno.gff3 ../`;
chdir "../";
}
#exit;

# make summary table for the annotation
my $genome_info = `perl $count_base $genome`;
my ($genome_size, $seq_count) = ((split /\s+/, $genome_info)[1] - (split /\s+/, $genome_info)[2], (split /\s+/, $genome_info)[4]);
`perl $gff2bed $genome.EDTA.TEanno.gff3 structural > $genome.EDTA.TEanno.bed`;
`perl -nle 'my (\$chr, \$s, \$e, \$anno, \$dir, \$supfam)=(split)[0,1,2,3,8,12]; print "10000 0.001 0.001 0.001 \$chr \$s \$e NA \$dir \$anno \$supfam"' $genome.EDTA.TEanno.bed > $genome.EDTA.TEanno.out`;

#`perl -nle 'next if /^#/; my (\$chr, undef, \$supfam, \$s, \$e, undef, \$dir, undef, \$anno)=(split); next if \$supfam=~/target_site_duplication|long_terminal_repeat/i; \$anno=\$1 if \$anno=~/Name=(.*?);/; print "10000 0.001 0.001 0.001 \$chr \$s \$e NA \$dir \$anno \$supfam"' $genome.EDTA.TEanno.gff3 > $genome.EDTA.TEanno.out`;

print "$genome_size\t$seq_count\n";
`perl $buildSummary -maxDiv 40 -genome_size $genome_size -seq_count $seq_count $genome.EDTA.TEanno.out > $genome.EDTA.TEanno.sum 2>/dev/null`;
