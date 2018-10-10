#!/bin/bash
export TEJAAS="${HOME}/tejaas/bin/tejaas"

BASE_DIR="/cbscratch/franco"
DATASET="gtex"
INPUTDIR="${BASE_DIR}/datasets/gtex"
OUTDIRPATH="/cbscratch/franco/tejaas_output/output2/"

export GTFILE="${INPUTDIR}/GTEx_450Indiv_genot_imput_info04_maf01_HWEp1E6_dbSNP135IDs_donorIDs_dosage_chr__CHROM__.gz"
export DONORS="${INPUTDIR}/donor_ids.fam"
export GENINF="${INPUTDIR}/gencode.v19.annotation.gtf.gz"

method="jpa-rr"

##################################################
################ WHOLE BLOOD #####################
##################################################
export GXFILE="${INPUTDIR}/Whole_Blood_Analysis.v6p.normalized.expression.txt"

# ~~ with cismask ~~
export EXTRAFLAGS="--cismask" #\t--force-cis # \t--force-trans
export OUTDIRBASE="${OUTDIRPATH}/${DATASET}_wb_cismask"

while read gtex_beta; do
  ./job_submit.sh gtex_master.jobsub $gtex_beta $method
done < sbetas.txt


# ~~ without cismask ~~
export EXTRAFLAGS=""
export OUTDIRBASE="${OUTDIRPATH}/${DATASET}_wb"

while read gtex_beta; do
  ./job_submit.sh gtex_master.jobsub $gtex_beta $method
done < sbetas.txt


##################################################
########### WHOLE BLOOD - NO PEER ################
##################################################
export GXFILE="${HOME}/tejaas/GTEx_wholeBlood_Normalzed_NoPEER_lmcorrected.txt"

# ~~ with cismask ~~
export EXTRAFLAGS="--cismask" #\t--force-cis # \t--force-trans
export OUTDIRBASE="${OUTDIRPATH}/${DATASET}_wb_cismask_nopeer"

while read gtex_beta; do
  ./job_submit.sh gtex_master.jobsub $gtex_beta $method
done < sbetas.txt


# ~~ without cismask ~~
export EXTRAFLAGS=""
export OUTDIRBASE="${OUTDIRPATH}/${DATASET}_wb_nopeer"

while read gtex_beta; do
  ./job_submit.sh gtex_master.jobsub $gtex_beta $method
done < sbetas.txt