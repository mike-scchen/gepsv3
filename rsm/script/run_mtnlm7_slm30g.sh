#!/bin/bash
#PJM --name "run_mtnlm7_slm30g"
#PJM --rsc-list "elapse=01:00:00"
#PJM --rsc-list "node=1"
#PJM --rsc-list "node-mem=24Gi"
#PJM --rsc-list "rscgrp=CWB"
#PJM --no-stging
#PJM -j

# Input files needed:
# fort.15  -> thirty.second.antarctic.new.bin
# fort.235 -> gmted2010.30sec.fine
#          -> gtopo30_gg.fine
# rsm-domain.nml

# Preset: 12km domain
# ARGUMENT="324 216"

# Preset:  5km domain
# ARGUMENT="768 432"

# Preset: Test domain
ARGUMENT="192 192"

# Custom:
# IGRD: Domain zonal dimension
# JGRD: Domain meridional dimension
# ARGUMENT="IGRD JGRD"

if [[ -z "$ARGUMENT" ]]; then
    echo "run_mtnlm7_slm30g: Must provide the 'ARGUMENT' environment variable!"
    exit
fi

rm -fv rsm-*.bin

echo "run_mtnlm7_slm30g: Using 'ARGUMENT=$ARGUMENT' as input"
./mtnlm7_slm30g.exe -Wl,-T15,-T235 <<< "$ARGUMENT"

