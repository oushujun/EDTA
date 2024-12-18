#!/usr/bin/env bash

set -euo pipefail

config_version=$(sed -n "/^\s*version\s*=\s*'/s/version//p" nextflow.config | tr -d "=[:space:]'")
perl_version=$(sed -n 's|^my $version = "\(.*\)";|\1|p' EDTA.pl | tr -d '[:space:]')

if [[ "v$config_version" != $perl_version ]]; then
    echo "config_version (v$config_version) != perl_version ($perl_version)"
    exit 1
fi
