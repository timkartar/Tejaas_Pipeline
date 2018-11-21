#!/bin/bash

source "./CONFIG_GNETLMM.sh"

###############################################################
# Compute anchors 
# We compute which genes have a genetic anchor
# ---> Here tune cis or trans analysis <---
# anchor_thresh is the p-value threshold used for calling an anchor association. 
# We recommend using a stringent threshold, i. e. genome-wide significance.
JOBNAME="anchors"
bsub -n 8 -a openmp -q mpi -W ${RUNTIME} -R span[hosts=1] -R cbscratch \
                    -J ${JOBNAME} \
                    -o ${OUTPATH}/${JOBNAME}.log \
                    -e ${OUTPATH}/${JOBNAME}.err \
                    $PYTHON $GNET_PATH/GNetLMM/bin/gNetLMM_analyse --compute_anchors  \
                    --bfile $BFILE --pfile $PFILE --assoc0file $ASSOC0FILE --anchorfile $ANCHORFILE \
                    --anchor_thresh=$ANCHOR_THRESH  --window=$WINDOW --cis
