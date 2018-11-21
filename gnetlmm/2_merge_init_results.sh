#!/bin/bash

source "./CONFIG_GNETLMM.sh"

# Merging results
$PYTHON $GNET_PATH/GNetLMM/bin/gNetLMM_analyse --merge_assoc0_scan  --assoc0file $ASSOC0FILE --nSnps $STEP --bfile $BFILE

# Merging results
$PYTHON $GNET_PATH/GNetLMM/bin/gNetLMM_analyse --merge_corr  --gfile $GFILE  --pfile $PFILE --nTraits $GSTEP
