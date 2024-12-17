#!/usr/bin/env bash

script_path=$(dirname $0)

cat > final/final_filter.pl << EOF

my \$genome             = 'genome';
my \$threads            = $1;
my \$rm_threads         = $1;

my \$overwrite          = 1;

my \$cleanup_nested     = "$script_path/cleanup_nested.pl";
my \$rename_TE          = "$script_path/rename_TE.pl";
my \$reclassify         = "$script_path/classify_by_lib_RM.pl";
my \$output_by_list     = '$script_path/output_by_list.pl';
my \$rename_by_list     = "$script_path/rename_by_list.pl";
my \$filter_gff         = "$script_path/filter_gff3.pl";
my \$format_gff3        = "$script_path/format_gff3.pl";
my \$add_id             = "$script_path/add_id.pl";

my \$blastplus          = ''; # Assumed on PATH
my \$repeatmasker       = ''; # Assumed on PATH

EOF

sed -n 607,612p "$script_path/EDTA.pl" \
    >> final/final_filter.pl

sed -n 629,676p "$script_path/EDTA.pl" \
    >> final/final_filter.pl
