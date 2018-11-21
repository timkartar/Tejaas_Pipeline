#!/bin/bash

source "./CONFIG_GNETLMM.sh"

###############################################################
# look at the number of Traits to analyze
# # Find v-structures
# corrThresh is the q-value threshold to call a gene-gene correlation significant (default: 0.01)
# indThresh is the q-value threshold to call a gene-gene correlation not significant (default: 0.1)
# max_genes can be used to set the maximum number of genes in the conditioning set (default: inf)
# nTraits is the number of traits to be analyzed (default: all)
# ¿¿?? if the SNP lies between gene_start-window and gene_end+window, the gene cannot be in the conditioning set of the marker.
for i in $(seq 0 $GSTEP $NTRAITS)
do
    JOBNAME="v_structures_${i}"
    bsub -n 8 -a openmp -q mpi -W ${RUNTIME} -R span[hosts=1] -R cbscratch \
                    -J ${JOBNAME} \
                    -o ${OUTPATH}/${JOBNAME}.log \
                    -e ${OUTPATH}/${JOBNAME}.err \
    $PYTHON $GNET_PATH/GNetLMM/bin/gNetLMM_analyse --find_vstructures  --pfile $PFILE  --gfile $GFILE --anchorfile $ANCHORFILE  \
    --assoc0file $ASSOC0FILE  --window $WINDOW --vfile $VFILE --bfile $BFILE --startTraitIdx $i --nTraits $GSTEP
done