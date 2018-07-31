#!/bin/bash

TEJAAS="${HOME}/tejaas/bin/tejaas"

CARDIO_INPUTDIR="${HOME}/datasets/cardiogenics"
CARDIO_OUTDIRBASE="${HOME}/pipeline/cardio_results"

CARDIO_GTFILE="${CARDIO_INPUTDIR}/genotypes/CG___CHROM__.imputed.gz"
CARDIO_GXFILE="${CARDIO_INPUTDIR}/cardio_mono_expr.txt"
CARDIO_DONORS="${CARDIO_INPUTDIR}/genotypes/CG.sample"
CARDIO_GENINF="${CARDIO_INPUTDIR}/gencode.v19.annotation.gtf.gz"




GTEX_INPUTDIR="${HOME}/datasets/gtex"
GTEX_OUTDIRBASE="${HOME}/pipeline/gtex_results"

GTEX_GTFILE="${GTEX_INPUTDIR}/GTEx_450Indiv_genot_imput_info04_maf01_HWEp1E6_dbSNP135IDs_donorIDs_dosage_chr__CHROM__.gz"
GTEX_GXFILE="${GTEX_INPUTDIR}/Whole_Blood_Analysis.v6p.normalized.expression.txt"
GTEX_DONORS="${GTEX_INPUTDIR}/donor_ids.fam"
GTEX_GENINF="${GTEX_INPUTDIR}/gencode.v19.annotation.gtf.gz"
