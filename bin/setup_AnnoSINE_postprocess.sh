#!/usr/bin/env bash

script_path=$(dirname $0)

cat > input/AnnoSINE_postprocess.pl << EOF
my \$genome             = 'genome';
my \$threads            = $1;

my \$cleanup_misclas    = '$script_path/cleanup_misclas.pl';
my \$cleanup_tandem     = '$script_path/cleanup_tandem.pl';

my \$TEsorter           = ''; # Assumed on PATH
my \$trf                = '$(which trf)';

EOF

sed -n 490,497p "$script_path/EDTA_raw.pl" \
    >> input/AnnoSINE_postprocess.pl
