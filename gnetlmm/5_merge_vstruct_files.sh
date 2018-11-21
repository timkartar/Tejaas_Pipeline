#!/bin/bash

source "./CONFIG_GNETLMM.sh"

# Merging csv files
$PYTHON $GNET_PATH/GNetLMM/bin/gNetLMM_postprocess --concatenate --infiles $VFILE      --outfile $VFILE
