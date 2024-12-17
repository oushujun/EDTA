#!/usr/bin/env bash

script_path=$(dirname $0)

cat > raw/combine_intact_TEs.pl << EOF
my \$genome             = 'genome';
my \$threads            = $1;

my \$bed2gff            = "$script_path/bed2gff.pl";

EOF

sed -n 495,504p "$script_path/EDTA.pl" \
    >> raw/combine_intact_TEs.pl