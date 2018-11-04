#!/bin/bash
# export TEJAAS="${HOME}/tejaas_stable/tejaas/bin/tejaas"
export TEJAAS="${HOME}/tejaas/bin/tejaas"

BASE_DIR="/cbscratch/franco"
DATASET="cardio"
INPUTDIR="${BASE_DIR}/datasets/cardiogenics"

# OUTDIRPATH="/cbscratch/franco/tejaas_output/output3"
OUTDIRPATH="/cbscratch/franco/tejaas_output/output3_it1000G"

export GTPATH="${INPUTDIR}/genotypes/CG___CHROM__.imputed.gz"
export DONORS="${INPUTDIR}/genotypes/CG.sample"
export GENINF="${INPUTDIR}/gencode.v19.annotation.gtf.gz"

method="jpa-rr"

run_cismask=true
run_nocismask=false
run_cismask_us=false

##################################################
################ MONOCYTES #######################
##################################################
# export GXFILE="${INPUTDIR}/cardio_mono_expr.txt"
export GXFILE="${INPUTDIR}/cardio_mono_matrix_nodup.txt.gencode_filter"


# ~~ with cismask ~~
if [ $run_cismask = true ]; then
    echo "Calculating monocytes with cismask"
    export EXTRAFLAGS="--cismask" #\t--force-cis # \t--force-trans
    export OUTDIRBASE="${OUTDIRPATH}/${DATASET}_mono_cismask"

    while read cardio_beta; do
        ./job_submit.sh cardio_master.jobsub $cardio_beta $method
    done < sbetas.txt
fi


# ~~ with cismask + user_samples ~~
if [ $run_cismask_us = true ]; then
    echo "Calculating monocytes with cismask + usersamples"
    export EXTRAFLAGS="--cismask" #\t--force-cis # \t--force-trans
    USERSAMPL="${HOME}/Tejaas_Pipeline/devtools/mono_usersamples.txt"
    export OPT_SAMPLES="--selected-donors\t${USERSAMPL}"
    export OUTDIRBASE="${OUTDIRPATH}/${DATASET}_mono_cismask_us"
    while read cardio_beta; do
        ./job_submit.sh cardio_master.jobsub $cardio_beta $method
    done < sbetas.txt
    export OPT_SAMPLES=""
fi


# ~~ without cismask ~~
if [ $run_nocismask = true ]; then
    echo "Calculating monocytes without cismask"
    export EXTRAFLAGS=""
    export OUTDIRBASE="${OUTDIRPATH}/${DATASET}_mono"
    while read cardio_beta; do
        ./job_submit.sh cardio_master.jobsub $cardio_beta $method
    done < sbetas.txt
fi


## For other gene expression sets, just change GXFILE
####################################################
################ MACROPHAGES #######################
####################################################
# export GXFILE="${INPUTDIR}/cardio_macro_expr.txt"
export GXFILE="${INPUTDIR}/cardio_macro_matrix_nodup.txt.gencode_filter"

# ~~ with cismask ~~
if [ $run_cismask = true ]; then
    echo "Calculating macrophages with cismask"
    export EXTRAFLAGS="--cismask"
    export OUTDIRBASE="${OUTDIRPATH}/${DATASET}_macro_cismask"
    while read cardio_beta; do
        ./job_submit.sh cardio_master.jobsub $cardio_beta $method
    done < sbetas.txt
fi


# ~~ with cismask + user_samples ~~
if [ $run_cismask_us = true ]; then
    echo "Calculating macrophages with cismask +  usersamples"
    export EXTRAFLAGS="--cismask"
    USERSAMPL="${HOME}/Tejaas_Pipeline/devtools/macro_usersamples.txt"
    export OPT_SAMPLES="--selected-donors\t${USERSAMPL}"
    export OUTDIRBASE="${OUTDIRPATH}/${DATASET}_macro_cismask_us"
    while read cardio_beta; do
        ./job_submit.sh cardio_master.jobsub $cardio_beta $method
    done < sbetas.txt
    export OPT_SAMPLES=""
fi

# ~~ without cismask ~~
if [ $run_nocismask = true ]; then
    echo "Calculating macrophages without cismask"
    export EXTRAFLAGS=""
    export OUTDIRBASE="${OUTDIRPATH}/${DATASET}_macro"
    while read cardio_beta; do
        ./job_submit.sh cardio_master.jobsub $cardio_beta $method
    done < sbetas.txt
fi

# GXFILE="${INPUTDIR}/cardio_mono_expr.txt.common_samples"
# GXFILE="${INPUTDIR}/cardio_macro_expr.txt.common_samples"
