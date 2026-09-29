#!/bin/bash
# Command-line access to the TOPMed Imputation Server.
#
# The TOPMed reference panel (r3, 133,597 samples) cannot be downloaded. It can
# only be used through the server (https://imputation.biodatacatalyst.nhlbi.nih.gov).
# This script installs imputationbot (github.com/lukfor/imputationbot), which
# submits jobs to the server and downloads the results through its REST API, so
# no browser or scp is needed.
#
# One-time manual step after running this script: create an API token (server ->
# profile -> API token), then run
#   bin/imputationbot add-instance
#   URL:   https://imputation.biodatacatalyst.nhlbi.nih.gov
#   token: <paste>
# imputationbot keeps the token in ~/.imputationbot (a small config file in home).
#
# For local Minimac4 runs, the panel is the NYGC 1000G 30x one, already on disk:
#   /u/project/cluo/terencew/demux_benchmark/pool_design/vcf/1000G/by_chrom/
# Eagle hg38 genetic map: /u/project/cluo/terencew/programs/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz
# minimac4 v4.1.6:        /u/project/cluo/terencew/programs/build/Minimac4/build/minimac4
#
# Usage: bash 00_setup_imputationbot.sh

set -euo pipefail

TOPMED=/u/project/cluo/terencew/reference/topmed
VERSION=2.1.0
ZIP=imputationbot-${VERSION}-linux.zip
URL=https://github.com/lukfor/imputationbot/releases/download/v${VERSION}/${ZIP}

mkdir -p "$TOPMED/bin"
cd "$TOPMED/bin"

if [ -x imputationbot ]; then
    echo "$(date): imputationbot already installed, skipping download"
else
    echo "$(date): downloading $URL"
    time curl -fL "$URL" -o "$ZIP"
    unzip -o "$ZIP"
    chmod +x imputationbot
fi

./imputationbot version
echo "$(date): done. Next: $TOPMED/bin/imputationbot add-instance (see header)"
