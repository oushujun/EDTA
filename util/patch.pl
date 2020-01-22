#!/usr/bin/env perl -w
use strict;
use FindBin;
use File::Basename;

my $usage = "\n\tPatch EDTA results from v1.6.x to v1.7.x. After this patch, rerun EDTA.pl -step final
		perl patch.pl genome.fa\n";

## preset
my $script_path = $FindBin::Bin;
my $rename_LTR = "$script_path/rename_LTR.pl";
my $output_by_list = "$script_path/output_by_list.pl";
my $make_gff = "$script_path/make_gff_with_intact.pl";
my $bed2gff = "$script_path/bed2gff.pl";
my $gff2bed = "$script_path/gff2bed.pl";

my $genome = $ARGV[0];

## raw/LTR
chdir "$genome.EDTA.raw/LTR";

# generate annotated output and gff
`ln -s ../../$genome`;
`perl $rename_LTR $genome $genome.LTR.intact.fa $genome.defalse > $genome.LTR.intact.fa.anno`;
`mv $genome.LTR.intact.fa.anno $genome.LTR.intact.fa` if -s "$genome.LTR.intact.fa.anno";
`cp $genome.LTR.intact.fa $genome.LTR.intact.fa.gff3 ../`;

## raw/TIR
chdir "../TIR";

# get gff of intact TIR elements
`perl -nle 's/\\-\\+\\-/_Len:/; my (\$chr, \$method, \$supfam, \$s, \$e, \$anno) = (split)[0,1,2,3,4,8]; my \$class='DNA'; \$class='MITE' if \$e-\$s+1 <= 600; \$anno = "ID=\$chr:\$s..\$e#\$class/\$supfam;\$anno"; print "\$chr \$method \$class/\$supfam \$s \$e . . . \$anno \$chr:\$s..\$e"' ./TIR-Learner-Result/TIR-Learner_FinalAnn.gff3| perl $output_by_list 10 - 1 $genome.TIR.raw.fa -MSU0 -MSU1 | awk '{\$10=""; print \$0}' | perl -nle 's/\\s+/\\t/g; print \$_' > $genome.TIR.intact.fa.gff`;
`perl -i -nle 's/_([0-9]+)_([0-9]+)#/:\$1..\$2#/; print \$_' $genome.TIR.intact.fa`;
`cp $genome.TIR.intact.fa.gff $genome.TIR.intact.fa ../`;

## raw/Helitron
chdir "../Helitron";
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

