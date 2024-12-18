#!/usr/bin/env bash

script_path=$(dirname $0)

cat > input/HelitronScanner_postprocess.pl << EOF
my \$genome             = 'genome';
my \$threads            = $1;

my \$flank_filter       = "$script_path/flanking_filter.pl";
my \$output_by_list     = '$script_path/output_by_list.pl';
my \$cleanup_tandem     = '$script_path/cleanup_tandem.pl';
my \$cleanup_misclas    = '$script_path/cleanup_misclas.pl';
my \$make_bed           = "$script_path/make_bed_with_intact.pl";
my \$bed2gff            = "$script_path/bed2gff.pl";

my \$blastplus          = ''; # Assumed on PATH
my \$mdust              = ''; # Assumed on PATH
my \$trf                = '$(which trf)';
my \$TEsorter           = ''; # Assumed on PATH

EOF

sed -n 705,722p "$script_path/EDTA_raw.pl" \
    >> input/HelitronScanner_postprocess.pl

sed -n 726,729p "$script_path/EDTA_raw.pl" \
    >> input/HelitronScanner_postprocess.pl
