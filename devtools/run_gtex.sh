#!/bin/bash
export TEJAAS="${HOME}/tejaas_stable/tejaas/bin/tejaas"

BASE_DIR="/cbscratch/franco"
DATASET="gtex"
INPUTDIR="${BASE_DIR}/datasets/gtex"
OUTDIRPATH="/cbscratch/franco/tejaas_output/output3/"

export GTPATH="${INPUTDIR}/GTEx_450Indiv_genot_imput_info04_maf01_HWEp1E6_dbSNP135IDs_donorIDs_dosage_chr__CHROM__.gz"
export DONORS="${INPUTDIR}/donor_ids.fam"
export GENINF="${INPUTDIR}/gencode.v19.annotation.gtf.gz"

method="jpa-rr"

run_cismask=true
run_nocismask=true
run_cismask_nopeer=false
run_nocismask_nopeer=false

wb=true
hlv=false
ms=false

##################################################
################ WHOLE BLOOD #####################
##################################################
if [ $wb = true ]; then
    export GXFILE="${INPUTDIR}/Whole_Blood_Analysis.v6p.normalized.expression.txt"
    TISSUE="wb"
    # ~~ with cismask ~~
    if [ $run_cismask = true ]; then
        echo "Calculating GTEx Whole Blood with cismask"
        export EXTRAFLAGS="--cismask" #\t--force-cis # \t--force-trans
        export OUTDIRBASE="${OUTDIRPATH}/${DATASET}_${TISSUE}_cismask"

        while read gtex_beta; do
          ./job_submit.sh gtex_master.jobsub $gtex_beta $method
        done < sbetas.txt
    fi


    # ~~ without cismask ~~
    if [ $run_nocismask = true ]; then
        echo "Calculating GTEx Whole Blood without cismask"
        export EXTRAFLAGS=""
        export OUTDIRBASE="${OUTDIRPATH}/${DATASET}_${TISSUE}"

        while read gtex_beta; do
          ./job_submit.sh gtex_master.jobsub $gtex_beta $method
        done < sbetas.txt
    fi


    ##################################################
    ########### WHOLE BLOOD - NO PEER ################
    ##################################################
    export GXFILE="${HOME}/tejaas/GTEx_wholeBlood_Normalzed_NoPEER_lmcorrected.txt"
    TISSUE="wb"
    # ~~ with cismask ~~
    if [ $run_cismask_nopeer = true ]; then
        echo "Calculating GTEx Whole Blood with cismask and no PEER"
        export EXTRAFLAGS="--cismask" #\t--force-cis # \t--force-trans
        export OUTDIRBASE="${OUTDIRPATH}/${DATASET}_${TISSUE}_cismask_nopeer"

        while read gtex_beta; do
          ./job_submit.sh gtex_master.jobsub $gtex_beta $method
        done < sbetas.txt
    fi


    # ~~ without cismask ~~
    if [ $run_nocismask_nopeer = true ]; then
        echo "Calculating GTEx Whole Blood without cismask and no PEER"
        export EXTRAFLAGS=""
        export OUTDIRBASE="${OUTDIRPATH}/${DATASET}_${TISSUE}_nopeer"

        while read gtex_beta; do
          ./job_submit.sh gtex_master.jobsub $gtex_beta $method
        done < sbetas.txt
    fi
fi


##################################################
################ Heart Left Ventricle ############
##################################################
if [ $hlv = true ]; then
    GTEx_expr_path="expression/GTEx_Analysis_v6p_eQTL_expression_matrices"
    export GXFILE="${INPUTDIR}/${GTEx_expr_path}/Heart_Left_Ventricle_Analysis.v6p.normalized.expression.bed.tejaas.txt.gencode_filter"
    TISSUE="hlv"
    # ~~ with cismask ~~
    if [ $run_cismask = true ]; then
        echo "Calculating GTEx Heart Left Ventricle with cismask"
        export EXTRAFLAGS="--cismask" #\t--force-cis # \t--force-trans
        export OUTDIRBASE="${OUTDIRPATH}/${DATASET}_${TISSUE}_cismask"

        while read gtex_beta; do
          ./job_submit.sh gtex_master.jobsub $gtex_beta $method
        done < sbetas.txt
    fi


    # ~~ without cismask ~~
    if [ $run_nocismask = true ]; then
        echo "Calculating GTEx Heart Left Ventricle without cismask"
        export EXTRAFLAGS=""
        export OUTDIRBASE="${OUTDIRPATH}/${DATASET}_${TISSUE}"

        while read gtex_beta; do
          ./job_submit.sh gtex_master.jobsub $gtex_beta $method
        done < sbetas.txt
    fi
fi


##################################################
################ Muscle Skeletal #################
##################################################
if [ $ms = true ]; then
    GTEx_expr_path="expression/GTEx_Analysis_v6p_eQTL_expression_matrices"
    export GXFILE="${INPUTDIR}/${GTEx_expr_path}/Muscle_Skeletal_Analysis.v6p.normalized.expression.bed.tejaas.txt.gencode_filter"
    TISSUE="ms"
    # ~~ with cismask ~~
    if [ $run_cismask = true ]; then
        echo "Calculating GTEx Muscle Skeletal with cismask"
        export EXTRAFLAGS="--cismask" #\t--force-cis # \t--force-trans
        export OUTDIRBASE="${OUTDIRPATH}/${DATASET}_${TISSUE}_cismask"

        while read gtex_beta; do
          ./job_submit.sh gtex_master.jobsub $gtex_beta $method
        done < sbetas.txt
    fi


    # ~~ without cismask ~~
    if [ $run_nocismask = true ]; then
        echo "Calculating GTEx Muscle Skeletal without cismask"
        export EXTRAFLAGS=""
        export OUTDIRBASE="${OUTDIRPATH}/${DATASET}_${TISSUE}"

        while read gtex_beta; do
          ./job_submit.sh gtex_master.jobsub $gtex_beta $method
        done < sbetas.txt
    fi
fi