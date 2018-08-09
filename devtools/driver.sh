#!/bin/bash
while read gtex_beta; do
	./gtex_job_submit.sh gtex_master.jobsub $gtex_beta
done < sbetas.txt

while read cardio_beta; do
    ./cardio_job_submit.sh cardio_master.jobsub $cardio_beta
done < sbetas.txt
