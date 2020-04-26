#!/usr/bin/env perl
use warnings;
use strict;
use FindBin;
use File::Basename;

my $usage = "\n\tPatch EDTA results from v1.6.x to v1.7.4+. After this patch, rerun EDTA.pl -step final
		perl patch.pl genome.fa\n";

## preset
my $script_path = $FindBin::Bin;
my $rename_LTR = "$script_path/rename_LTR.pl";
my $output_by_list = "$script_path/output_by_list.pl";
my $make_gff = "$script_path/make_gff_with_intact.pl";
my $bed2gff = "$script_path/bed2gff.pl";
my $gff2bed = "$script_path/gff2bed.pl";

my $genome = $ARGV[0];
my $genome_ori = $genome;

## raw/LTR
`cp $genome $genome.mod`;
`mv $genome.EDTA.raw $genome.mod.EDTA.raw 2>/dev/null`;
chdir "$genome.mod.EDTA.raw/LTR";
`for i in \`ls\`; do mv -f \$i \$(echo \$i | perl -nle 's/$genome_ori/$genome_ori.mod/g; s/.mod.mod/.mod/g; print \$_') 2>/dev/null; done`;
$genome = "$genome.mod";

# generate annotated output and gff
`ln -s ../../$genome` unless -s $genome;
`perl $rename_LTR $genome $genome.LTR.intact.fa $genome.defalse > $genome.LTR.intact.fa.anno`;
`mv $genome.LTR.intact.fa.anno $genome.LTR.intact.fa` if -s "$genome.LTR.intact.fa.anno";
`cp $genome.LTR.intact.fa $genome.LTR.intact.fa.gff3 ../`;

## raw/TIR
chdir "../TIR";
`for i in \`ls\`; do mv -f \$i \$(echo \$i | perl -nle 's/$genome_ori/$genome_ori.mod/g; s/.mod.mod/.mod/g; print \$_') 2>/dev/null; done`;

# get gff of intact TIR elements
`perl -nle 's/\\-\\+\\-/_Len:/; my (\$chr, \$method, \$supfam, \$s, \$e, \$anno) = (split)[0,1,2,3,4,8]; my \$class='DNA'; \$class='MITE' if \$e-\$s+1 <= 600; \$anno = "ID=\$chr:\$s..\$e#\$class/\$supfam;\$anno"; print "\$chr \$method \$class/\$supfam \$s \$e . . . \$anno \$chr:\$s..\$e"' ./TIR-Learner-Result/TIR-Learner_FinalAnn.gff3| perl $output_by_list 10 - 1 $genome.TIR.raw.fa -MSU0 -MSU1 | awk '{\$10=""; print \$0}' | perl -nle 's/\\s+/\\t/g; print \$_' > $genome.TIR.intact.fa.gff`;
`perl -nle 's/_([0-9]+)_([0-9]+)#/:\$1..\$2#/; print \$_' $genome.TIR.raw.fa > $genome.TIR.intact.fa`;
`cp $genome.TIR.intact.fa.gff $genome.TIR.intact.fa ../`;

## raw/Helitron
chdir "../Helitron";
`for i in \`ls\`; do mv -f \$i \$(echo \$i | perl -nle 's/$genome_ori/$genome_ori.mod/g; s/.mod.mod/.mod/g; print \$_') 2>/dev/null; done`;
`perl -nle 's/>(.*)/>\$1#DNA\\/Helitron/ unless /Helitron/; print \$_' $genome.Helitron.raw.fa > $genome.Helitron.intact.fa`;
`perl $make_gff $genome.Helitron.intact.fa`;
`cp $genome.Helitron.intact.fa $genome.Helitron.intact.fa.gff ../`;

## raw/
chdir '../';

# combine intact TEs
`cat $genome.LTR.intact.fa $genome.TIR.intact.fa $genome.Helitron.intact.fa > $genome.EDTA.intact.fa`;
`cat $genome.LTR.intact.fa.gff3 $genome.TIR.intact.fa.gff $genome.Helitron.intact.fa.gff | perl $gff2bed - structural > $genome.EDTA.intact.bed`;
`perl $bed2gff $genome.EDTA.intact.bed`;
`mv $genome.EDTA.intact.bed.gff $genome.EDTA.intact.gff`;
`cp $genome.EDTA.intact.gff ../`;

## rename filtered
chdir '../';
`mv $genome_ori.EDTA.combine $genome.EDTA.combine`;
chdir "$genome.EDTA.combine";
`for i in \`ls\`; do mv -f \$i \$(echo \$i | perl -nle 's/$genome_ori/$genome_ori.mod/g; s/.mod.mod/.mod/g; print \$_') 2>/dev/null; done`;

## rename final
chdir '../';
`mv $genome_ori.EDTA.final $genome.EDTA.final`;
chdir "$genome.EDTA.final";
`for i in \`ls\`; do mv -f \$i \$(echo \$i | perl -nle 's/$genome_ori/$genome_ori.mod/g; s/.mod.mod/.mod/g; print \$_') 2>/dev/null; done`;


