#!/bin/bash

TEJAAS="${HOME}/tejaas/bin/tejaas"
BASE_DIR="/cbscratch/franco"

 ##### Cardiogenics

DATASET="cardio"

# OUTDIRMOD="${DATASET}_forcecis_cismask_mpi"
OUTDIRNAME="${DATASET}_cismask_mpi_newtest"
# OUTDIRMOD="${DATASET}_macro_us" #user samples + cismask
# OUTDIRMOD="${DATASET}_mono_us" #user samples + cismask


INPUTDIR="${BASE_DIR}/datasets/cardiogenics"
OUTDIRBASE="/cbscratch/franco/tejaas_output/output/${OUTDIRNAME}"
GTFILE="${INPUTDIR}/genotypes/CG___CHROM__.imputed.gz"
DONORS="${INPUTDIR}/genotypes/CG.sample"
GENINF="${INPUTDIR}/gencode.v19.annotation.gtf.gz"

GXFILE="${INPUTDIR}/cardio_mono_expr.txt"
# GXFILE="${INPUTDIR}/cardio_macro_expr.txt"
# GXFILE="${INPUTDIR}/cardio_mono_expr.txt.common_samples"
# GXFILE="${INPUTDIR}/cardio_macro_expr.txt.common_samples"





DATASET="gtex"
OUTDIRNAME="${DATASET}_cismask_mpi_newtest"
INPUTDIR="${BASE_DIR}/datasets/gtex"
OUTDIRBASE="/cbscratch/franco/tejaas_output/output/${OUTDIRNAME}"
GTFILE="${INPUTDIR}/GTEx_450Indiv_genot_imput_info04_maf01_HWEp1E6_dbSNP135IDs_donorIDs_dosage_chr__CHROM__.gz"
GXFILE="${INPUTDIR}/Whole_Blood_Analysis.v6p.normalized.expression.txt"
DONORS="${INPUTDIR}/donor_ids.fam"
GENINF="${INPUTDIR}/gencode.v19.annotation.gtf.gz"
