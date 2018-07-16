#!/bin/bash
while read gtex_beta; do
	sed "s/_BETA_/${gtex_beta}/g" gtex_master.jobsub > gtex.jobsub
	./gtex_job_submit.sh gtex.jobsub $gtex_beta
done < sbetas.txt

while read cardio_beta; do
        sed "s/_BETA_/${cardio_beta}/g" cardio_master.jobsub > cardio.jobsub
        ./cardio_job_submit.sh cardio.jobsub $cardio_beta
done < sbetas.txt
