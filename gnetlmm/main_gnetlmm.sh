#!/bin/bash 

PYTHON3="/usr/users/fsimone/opt/miniconda/3/envs/env3.6/bin/python"
PYTHON_DIR="/usr/users/fsimone/opt/miniconda/3/envs/env2.7"
PYTHON="${PYTHON_DIR}/bin/python"
OUTPATH="/cbscratch/franco/tejaas_output/gnetlmm_112018_ms"
# Create outpath
if [ ! -d "$OUTPATH" ]; then
    mkdir -p "$OUTPATH"
    mkdir -p "$OUTPATH/out"
fi

LOCALDIR="/usr/users/fsimone/Tejaas_Pipeline/gnetlmm"

GNET_PATH="${HOME}/GNetLMM"
# BFILE="${GNET_PATH}/data/1000G_chr22/chrom22_subsample20_maf0.10" #specify here bed basename
BFILE="/cbscratch/franco/datasets/gtex/genotypes/vcfs_split_ms/bed/GTEx_nomissing_chr22" #specify here bed basename

# phenotype matrix file. Expression values are saved in pfile.matrix. 
# Each row corresponds to one gene, each column to one sample. 
# Row informations are saved in pfile.rows and column informations in pfile.cols. 
# Row and column files must have a header describing the information, while the expression matrix does not have a header. 
# The row file must contain the information gene_start, gene_end, gene_chrom and gene_ids. 
# Note that the first time the phenotype file is loaded, an index file pfile.cache.npy is created allowing faster access later on. 
# The samples in pfile and bfile must have the same ordering in order for the software to work.


# Generate phenotypes
# PFILE="$LOCALDIR/out/pheno"
# $PYTHON $GNET_PATH/GNetLMM/bin/gNetLMM_simPheno --bfile $BFILE --pfile $PFILE

PFILE="$LOCALDIR/expr_gencode_filtered/ms_gtex.35ncov_PEER_residuals.txt"

NLINES=`cat $PFILE.cols | wc -l`
NDONORS=$(echo $(( NLINES - 1)))

# ffile is the filename of the covariates file.
# FFILE="${GNET_PATH}/data/1000G_chr22/ones.txt"
FFILE="$LOCALDIR/ms.covariates.dummy.gnet"
$PYTHON3 format_covariates_gnetlmm.py --make-dummy --donors $NDONORS --output $FFILE


# cfile is the filename of the (genotype!) covariance matrix file.
CFILE="${OUTPATH}/covmat/chrom22"

# assoc0file is the base filename to which the results are written out. 
# The file assoc0file.pv.matrix contains the p-values, the file assoc0file.beta.matrix the SNP weights. 
# Each row corresponds to one marker, each column to one gene.



GFILE="${OUTPATH}/out/genes"


# we compute which genes have a genetic anchor
# anchorfile is the results filename (tab-delimited format) having the following columns: 
#   gene_ids, snp_ids, gene index, snp index and p-value.

ANCHOR_THRESH=1e-6
WINDOW=2000
ANCHORFILE="${OUTPATH}/out/cisanchor_thresh${ANCHOR_THRESH}_wnd${WINDOW}.txt"
VFILE="${OUTPATH}/out/vstructures_thresh${ANCHOR_THRESH}_wnd${WINDOW}"
ASSOCFILE="${OUTPATH}/out/gnetlmm_thresh${ANCHOR_THRESH}_wnd${WINDOW}"

# add paths for gNetLMM to the environment
echo $GNET_PATH > "${PYTHON_DIR}/lib/python2.7/site-packages/gnetlmm.pth"



PLOTFILE="${OUTPATH}/out/power.pdf"



# Compute covariance matrix
$GNET_PATH/GNetLMM/bin/gNetLMM_preprocess --compute_covariance --bfile $BFILE --cfile $CFILE


###############################################################
# Run initial association scan
# adjust number of SNPs to run scan to chromosome being used (right now is 40000 for testing)

# This took a lot! Paralelize!!
# Running initial scan
# ... finished in 4256.26197481 seconds

RUNTIME="6:00"
ASSOC0FILE="${OUTPATH}/out/lmm"
NSNPS_CHR=`cat $BFILE.bim |wc -l`
echo $NSNPS_CHR
STEP=2000
for i in $(seq 0 $STEP $NSNPS_CHR)
do
    JOBNAME="init_assoc_chrom22_${i}"

    echo bsub -n 8 -a openmp -q mpi -W ${RUNTIME} -R span[hosts=1] -R cbscratch \
                    -J ${JOBNAME} \
                    -o ${OUTPATH}/${JOBNAME}.log \
                    -e ${OUTPATH}/${JOBNAME}.err \
                    $PYTHON $GNET_PATH/GNetLMM/bin/gNetLMM_analyse --initial_scan --bfile $BFILE \
                                            --pfile $PFILE --cfile $CFILE.cov \
                                            --assoc0file $ASSOC0FILE --startSnpIdx $i --nSnps $STEP 
                                            ## --ffile $FFILE ## no covariates file, gx is already corrected
done

# /usr/users/fsimone/opt/miniconda/3/envs/env2.7/bin/python /usr/users/fsimone/GNetLMM/GNetLMM/bin/gNetLMM_analyse --initial_scan --bfile /cbscratch/franco/datasets/gtex/genotypes/vcfs_split_ms/bed/GTEx_chr22 --pfile /usr/users/fsimone/Tejaas_Pipeline/gnetlmm/expr_gencode_filtered/ms_gtex.35ncov_PEER_residuals.txt --cfile /cbscratch/franco/tejaas_output/gnetlmm_112018_ms/covmat/chrom22.cov --assoc0file /cbscratch/franco/tejaas_output/gnetlmm_112018_ms/out/lmm2 --startSnpIdx 0 --nSnps 2000 --ffile /usr/users/fsimone/Tejaas_Pipeline/gnetlmm/ms.covariates.dummy.gnet

# Merging results
$PYTHON $GNET_PATH/GNetLMM/bin/gNetLMM_analyse --merge_assoc0_scan  --assoc0file $ASSOC0FILE --nSnps $STEP --bfile $BFILE


# this can be run independently than the initial association analysis above
###############################################################
# Compute marginal gene-gene correlations
# nTraits is the number of traits analyzed in each partial file.
# pfile is the basename of the phenotype file. It is needed to get the total number of traits.
NTRAITS=`cat $PFILE.matrix |wc -l`
GSTEP=2000
for i in $(seq 0 $GSTEP $NTRAITS)
do
    JOBNAME="genegene_corr_${i}"
    bsub -n 8 -a openmp -q mpi -W ${RUNTIME} -R span[hosts=1] -R cbscratch \
                    -J ${JOBNAME} \
                    -o ${OUTPATH}/${JOBNAME}.log \
                    -e ${OUTPATH}/${JOBNAME}.err \
    $PYTHON $GNET_PATH/GNetLMM/bin/gNetLMM_analyse --gene_corr --pfile $PFILE --gfile $GFILE  --startTraitIdx $i --nTraits $GSTEP
done

# Merging results
$PYTHON $GNET_PATH/GNetLMM/bin/gNetLMM_analyse --merge_corr  --gfile $GFILE  --pfile $PFILE --nTraits $GSTEP



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
                    echo "$PYTHON $GNET_PATH/GNetLMM/bin/gNetLMM_analyse --compute_anchors  \
                    --bfile $BFILE --pfile $PFILE --assoc0file $ASSOC0FILE --anchorfile $ANCHORFILE \
                    --anchor_thresh=$ANCHOR_THRESH  --window=$WINDOW --cis"


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

# Merging csv files
$PYTHON $GNET_PATH/GNetLMM/bin/gNetLMM_postprocess --concatenate --infiles $VFILE      --outfile $VFILE


###############################################################
# update the associations for which we found vstructures:
# HEAVY STEP! # Update associationss
for i in $(seq 0 10 90)
do
     $PYTHON $GNET_PATH/GNetLMM/bin/gNetLMM_analyse --update_assoc --bfile $BFILE --pfile $PFILE --cfile $CFILE.cov --ffile $FFILE --vfile $VFILE --assocfile $ASSOCFILE --startTraitIdx $i --nTraits 10
done


# Merging csv files
$PYTHON $GNET_PATH/GNetLMM/bin/gNetLMM_postprocess --concatenate --infiles $ASSOCFILE  --outfile $ASSOCFILE

# Write to matrix
$PYTHON $GNET_PATH/GNetLMM/bin/gNetLMM_postprocess --merge_assoc --assoc0file $ASSOC0FILE --assocfile $ASSOCFILE



###############################################################
# Alternative to above!
# The arguments have the same specification as for --update_assoc, but the anchor gene is added as covariate to the model. 
# The association between the SNP and the focal gene is thereby blocked by the anchor gene. 
# If the snp - focal gene association is mediated by the anchor gene, conditioning on it will decrease the association strength.
# HEAVY STEP! # Block associationss
for i in $(seq 0 10 90)
do
     $PYTHON $GNET_PATH/GNetLMM/bin/gNetLMM_analyse --block_assoc --bfile $BFILE --pfile $PFILE --cfile $CFILE.cov --ffile $FFILE --vfile $VFILE --assocfile $ASSOCFILE.block --startTraitIdx $i --nTraits 10
done

# Merging csv files
$PYTHON $GNET_PATH/GNetLMM/bin/gNetLMM_postprocess --concatenate --infiles $ASSOCFILE.block  --outfile $ASSOCFILE.block

# Write to matrix
$PYTHON $GNET_PATH/GNetLMM/bin/gNetLMM_postprocess --merge_assoc --assoc0file $ASSOC0FILE --assocfile $ASSOCFILE.block



###############################################################
# Plot results
# This method does only work for simulated phenotypes as it plots the false positive rate vs. the true positive rate. 
# The true underlying gene-gene network must be deposited in pfile.Agene.
$PYTHON $GNET_PATH/GNetLMM/bin/gNetLMM_postprocess --plot_power --assocfile $ASSOCFILE --assoc0file $ASSOC0FILE --plotfile $PLOTFILE --pfile $PFILE --bfile $BFILE --window $WINDOW --blockfile $ASSOCFILE.block



###############################################################
# Creating nice output file for v-structures
# outfile is the basename of the nicely formatted vstructure results file. 
# Each updated association corresponds to one row. 
# The first column contains the anchor snp, 
# the second column the anchor gene, 
# the third column the focal gene, 
# the fourth column the orthogonal gene(s), 
# the fifth column the p-value without conditioning (vanilla LMM), 
# the sixth column with conditioning on the orthogonal genes (GNet-LMM), 
# the seventh column with conditioning on the orthogonal genes and the anchor gene (Block-LMM).
$PYTHON $GNET_PATH/GNetLMM/bin/gNetLMM_postprocess --nice_output --bfile $BFILE --pfile $PFILE --vfile $VFILE --assoc0file $ASSOC0FILE --assocfile $ASSOCFILE --blockfile $ASSOCFILE.block --outfile $ASSOCFILE.nice