#!/bin/bash

PYTHON3="/usr/users/fsimone/opt/miniconda/3/envs/env3.6/bin/python"
# gnetlmm requires python2.7
PYTHON_DIR="/usr/users/fsimone/opt/miniconda/3/envs/env2.7"
PYTHON="${PYTHON_DIR}/bin/python"

# current directory (for dev purposes)
LOCALDIR="/usr/users/fsimone/Tejaas_Pipeline/gnetlmm"

# where to write outputs
OUTPATH="/cbscratch/franco/tejaas_output/gnetlmm_112018_ms"



GNET_PATH="${HOME}/GNetLMM"
# bed genotype basename file (.bed .bim .fam)
BFILE="/cbscratch/franco/datasets/gtex/genotypes/vcfs_split_ms/bed/GTEx_nomissing_chr22" #specify here bed basename
PFILE="$LOCALDIR/expr_gencode_filtered/ms_gtex.35ncov_PEER_residuals.txt"

NLINES=`cat $PFILE.cols | wc -l`
NDONORS=$(echo $(( NLINES - 1)))

# ffile is the filename of the covariates file.
# we don't need them for now, the program automatically generates a vector of 1's
# FFILE="$LOCALDIR/ms.covariates.dummy.gnet"
# $PYTHON3 format_covariates_gnetlmm.py --make-dummy --donors $NDONORS --output $FFILE

# output filename of the genotype samples covariance matrix file (N x N matrix N= nº of samples).
CFILE="${OUTPATH}/covmat/chrom22"

# output filename for gene-gene correlations
GFILE="${OUTPATH}/out/genes"

# initial cis-scan options
ANCHOR_THRESH=1e-6

# window size of cis-snps, should be 1Mb but for dev a small window is fine for testing
WINDOW=2000
ANCHORFILE="${OUTPATH}/out/cisanchor_thresh${ANCHOR_THRESH}_wnd${WINDOW}.txt"
VFILE="${OUTPATH}/out/vstructures_thresh${ANCHOR_THRESH}_wnd${WINDOW}"
ASSOCFILE="${OUTPATH}/out/gnetlmm_thresh${ANCHOR_THRESH}_wnd${WINDOW}"
PLOTFILE="${OUTPATH}/out/power.pdf"

# add paths for gNetLMM to the environment
echo $GNET_PATH > "${PYTHON_DIR}/lib/python2.7/site-packages/gnetlmm.pth"

# init assoc analysis
RUNTIME="6:00"
ASSOC0FILE="${OUTPATH}/out/lmm"
NSNPS_CHR=`cat $BFILE.bim |wc -l`
STEP=2000

# gene-gene corr analysis
NTRAITS=`cat $PFILE.matrix |wc -l`
GSTEP=2000

