#!/bin/bash

source "./CONFIG_GNETLMM.sh"


###############################################################
# Run initial association scan
# adjust number of SNPs to run scan to chromosome being used (right now is 40000 for testing)

for i in $(seq 0 $STEP $NSNPS_CHR)
do
    JOBNAME="init_assoc_chrom22_${i}"

    bsub -n 8 -a openmp -q mpi -W ${RUNTIME} -R span[hosts=1] -R cbscratch \
                    -J ${JOBNAME} \
                    -o ${OUTPATH}/${JOBNAME}.log \
                    -e ${OUTPATH}/${JOBNAME}.err \
                    $PYTHON $GNET_PATH/GNetLMM/bin/gNetLMM_analyse --initial_scan --bfile $BFILE \
                                            --pfile $PFILE --cfile $CFILE.cov \
                                            --assoc0file $ASSOC0FILE --startSnpIdx $i --nSnps $STEP 
                                            ## --ffile $FFILE ## no covariates file, gx is already corrected
done



# this can be run independently than the initial association analysis above
###############################################################
# Compute marginal gene-gene correlations
# nTraits is the number of traits analyzed in each partial file.
# pfile is the basename of the phenotype file. It is needed to get the total number of traits.

for i in $(seq 0 $GSTEP $NTRAITS)
do
    JOBNAME="genegene_corr_${i}"
    bsub -n 8 -a openmp -q mpi -W ${RUNTIME} -R span[hosts=1] -R cbscratch \
                    -J ${JOBNAME} \
                    -o ${OUTPATH}/${JOBNAME}.log \
                    -e ${OUTPATH}/${JOBNAME}.err \
    $PYTHON $GNET_PATH/GNetLMM/bin/gNetLMM_analyse --gene_corr --pfile $PFILE --gfile $GFILE  --startTraitIdx $i --nTraits $GSTEP
done