#!/usr/bin/env bash

set -euo pipefail

EDTA_raw_pl_version=$(md5sum EDTA_raw.pl | cut -f1 -d' ')
bin_EDTA_raw_pl_version=$(md5sum bin/EDTA_raw.pl | cut -f1 -d' ')

if [[ $EDTA_raw_pl_version != $bin_EDTA_raw_pl_version ]]; then
    echo 'EDTA_raw.pl != bin/EDTA_raw.pl'
    echo 'Please synchronize bin/EDTA_raw.pl with EDTA_raw.pl'
    exit 1
fi
