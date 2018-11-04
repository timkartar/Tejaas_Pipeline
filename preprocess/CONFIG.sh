#!/bin/bash

TISSUEFILE="tissues.table.txt"
ENV=/home/fsimone/myenv/bin

# REQUIRED FILES!! download them from GTEx (for pheno file you need access)
RPKMFILE="GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz"
PHENOFILE="phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt"
READSFILE="GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct.gz"
DONORFILE="../donor_ids.fam"


RPKMOUTDIR="rpkms"              # tissue specific rpkms output dir
EXPROUTDIR="normalized_expr"    # normalized, uncorrected expr dir

GTEXCOVDIR="GTEx_Analysis_v6p_eQTL_covariates"  # dir for covariates downloaded from GTEx (public access)
COVDIR="PEER_covariates"                        # PEER covariates outdir
RUNPEER=false
NCOV=35                     # Nº of hidden confounders for PEER

LMOUTDIR="norm_lmcorrected" # outdir of linear model corrected expressions (PC+platform+sex, no peer!)

# If you want to correct for PEER, use _residuals output directly from step 2 script.