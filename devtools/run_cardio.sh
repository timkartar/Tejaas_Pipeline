#!/bin/bash
export TEJAAS="${HOME}/tejaas/bin/tejaas"

BASE_DIR="/cbscratch/franco"
DATASET="cardio"
INPUTDIR="${BASE_DIR}/datasets/cardiogenics"
OUTDIRPATH="/cbscratch/franco/tejaas_output/output/"

export GTFILE="${INPUTDIR}/genotypes/CG___CHROM__.imputed.gz"
export DONORS="${INPUTDIR}/genotypes/CG.sample"
export GENINF="${INPUTDIR}/gencode.v19.annotation.gtf.gz"

method="jpa-rr"

##################################################
################ MONOCYTES #######################
##################################################
export GXFILE="${INPUTDIR}/cardio_mono_expr.txt"

# ~~ with cismask ~~
export EXTRAFLAGS="--cismask" #\t--force-cis # \t--force-trans
export OUTDIRBASE="${OUTDIRPATH}/${DATASET}_mono_cismask"

while read cardio_beta; do
    ./job_submit.sh cardio_master.jobsub $cardio_beta $method
done < sbetas.txt

# ~~ with cismask + user_samples ~~

USERSAMPL="${HOME}/Tejaas_Pipeline/devtools/mono_usersamples.txt"
export OPT_SAMPLES="--selected-donors\t${USERSAMPL}"
export OUTDIRBASE="${OUTDIRPATH}/${DATASET}_mono_cismask_us"
while read cardio_beta; do
    ./job_submit.sh cardio_master.jobsub $cardio_beta $method
done < sbetas.txt

export OPT_SAMPLES=""


# ~~ without cismask ~~

export EXTRAFLAGS=""
export OUTDIRBASE="${OUTDIRPATH}/${DATASET}_mono"

while read cardio_beta; do
    ./job_submit.sh cardio_master.jobsub $cardio_beta $method
done < sbetas.txt


## For other gene expression sets, just change GXFILE
####################################################
################ MACROPHAGES #######################
####################################################
export GXFILE="${INPUTDIR}/cardio_macro_expr.txt"

# ~~ with cismask ~~

export EXTRAFLAGS="--cismask"
export OUTDIRBASE="${OUTDIRPATH}/${DATASET}_macro_cismask"

while read cardio_beta; do
    ./job_submit.sh cardio_master.jobsub $cardio_beta $method $mod
done < sbetas.txt


# ~~ with cismask + user_samples ~~

USERSAMPL="${HOME}/Tejaas_Pipeline/devtools/macro_usersamples.txt"
export OPT_SAMPLES="--selected-donors\t${USERSAMPL}"
export OUTDIRBASE="${OUTDIRPATH}/${DATASET}_macro_cismask_us"
while read cardio_beta; do
    ./job_submit.sh cardio_master.jobsub $cardio_beta $method
done < sbetas.txt
export OPT_SAMPLES=""


# ~~ without cismask ~~

export EXTRAFLAGS=""
export GXFILE="${INPUTDIR}/cardio_macro_expr.txt"
export OUTDIRBASE="${OUTDIRPATH}/${DATASET}_macro"

while read cardio_beta; do
    ./job_submit.sh cardio_master.jobsub $cardio_beta $method $mod
done < sbetas.txt

# GXFILE="${INPUTDIR}/cardio_mono_expr.txt.common_samples"
# GXFILE="${INPUTDIR}/cardio_macro_expr.txt.common_samples"