#!/bin/bash
# Create outpath
if [ ! -d "$OUTPATH" ]; then
    mkdir -p "$OUTPATH"
    mkdir -p "$OUTPATH/out"
fi

source "./CONFIG_GNETLMM.sh"

# Compute covariance matrix
$GNET_PATH/GNetLMM/bin/gNetLMM_preprocess --compute_covariance --bfile $BFILE --cfile $CFILE
