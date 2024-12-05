#!/usr/bin/env bash

script_path=$(dirname $0)

cat > input/LTR_retriever_postprocess.pl << EOF

my \$genome         = 'genome';
my \$threads        = $1;

my \$call_seq       = '$script_path/call_seq_by_list.pl';
my \$rename_LTR     = '$script_path/rename_LTR_skim.pl';
my \$mdust          = ''; # Assumed on PATH
my \$cleanup_tandem = '$script_path/cleanup_tandem.pl';
my \$trf            = '$(which trf)';
my \$TEsorter       = ''; # Assumed on PATH
my \$cleanup_misclas= '$script_path/cleanup_misclas.pl';
my \$output_by_list = '$script_path/output_by_list.pl';
my \$filter_gff     = '$script_path/filter_gff3.pl';

EOF

sed -n 443,462p "$script_path/EDTA_raw.pl" \
    >> input/LTR_retriever_postprocess.pl

sed -n 466,472p "$script_path/EDTA_raw.pl" \
    >> input/LTR_retriever_postprocess.pl
