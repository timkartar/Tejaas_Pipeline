#!/bin/bash

mkdir -p ../plots
while read gtex_beta; do
	while read cardio_beta; do
		python validate.py  $gtex_beta $cardio_beta
	done < sbetas.txt
done < sbetas.txt

