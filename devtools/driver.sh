#!/bin/bash

# USERSAMPL="${HOME}/Tejaas_Pipeline/devtools/macro_usersamples.txt"
USERSAMPL="${HOME}/Tejaas_Pipeline/devtools/mono_usersamples.txt"

# --force-cis # --force-trans

method="jpa-rr"
mod="--cismask" #\t--selected-donors\t${USERSAMPL}" #\t--force-cis"

# while read gtex_beta; do
# 	./gtex_job_submit.sh gtex_master.jobsub $gtex_beta $method $mod
# done < gtex_sbetas.txt

while read cardio_beta; do
    ./cardio_job_submit.sh cardio_master.jobsub $cardio_beta $method $mod
done < cardio_sbetas.txt
