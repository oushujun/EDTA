#!/usr/bin/env bash

script_path=$(dirname $0)

cat > input/RepeatModeler_postprocess.pl << EOF
my \$genome             = 'genome';
my \$threads            = $1;

my \$cleanup_misclas    = '$script_path/cleanup_misclas.pl';
my \$cleanup_tandem     = '$script_path/cleanup_tandem.pl';
my \$output_by_list     = '$script_path/output_by_list.pl';

my \$TEsorter           = ''; # Assumed on PATH
my \$trf                = '$(which trf)';

EOF

sed -n 567,579p "$script_path/EDTA_raw.pl" \
    >> input/RepeatModeler_postprocess.pl

sed -n 585,584p "$script_path/EDTA_raw.pl" \
    >> input/RepeatModeler_postprocess.pl
