#!/bin/bash

source paths.sh

master=$1
beta=$2

INPUTDIR="${GTEX_INPUTDIR}"
OUTDIRBASE="${GTEX_OUTDIRBASE}${beta}"
JOBSUBDIR="${OUTDIRBASE}/jobsub"
#TEJAAS="${HOME}/tejaas_mkl/bin/tejaas"
CWD=`pwd`

MAX_NSNP=20000
SNPCUT=0.001
GENCUT=0.001

if [ ! -d ${JOBSUBDIR} ]; then
    mkdir -p ${JOBSUBDIR}
fi
    
while read j; do

    echo "Submitting jobs for Chromosome ${j}."

    GTFILE="${INPUTDIR}/GTEx_450Indiv_genot_imput_info04_maf01_HWEp1E6_dbSNP135IDs_donorIDs_dosage_chr${j}.gz"
    GXFILE="${INPUTDIR}/Whole_Blood_Analysis.v6p.normalized.expression.txt"
    DONORS="${INPUTDIR}/donor_ids.fam"
    GENINF="${INPUTDIR}/gencode.v19.annotation.gtf.gz"
    JOBPREFIX="chr${j}"

    ## do not change below

    OUTDIR="${OUTDIRBASE}/chr${j}"
    
    if [ ! -d ${OUTDIR} ]; then
        mkdir -p ${OUTDIR}
    fi   

    TOTALSNPS=`zcat ${GTFILE} | wc -l`
    NJOBS=$(echo $(( TOTALSNPS/MAX_NSNP )))
    
    for (( i=0; i <= ${NJOBS}; i++ )); do 
        INDEX=`echo $i | awk '{printf "%03d", $1}'`
        JOBNAME="${JOBPREFIX}_${INDEX}"

        STARTSNP=$(( MAX_NSNP * i + 1 ))
        ENDSNP=$(( MAX_NSNP * (i+1) ))
        if [ $ENDSNP -gt $TOTALSNPS ]; then
            ENDSNP=${TOTALSNPS}
        fi
        INCSNP="${STARTSNP}:${ENDSNP}"

        OUTPRF="${OUTDIR}/chunk${INDEX}"

        # create the job submission file
        sed "s|_JOBNAME|${JOBNAME}|g;
             s|_GTFILE_|${GTFILE}|g;
             s|_FAMFIL_|${DONORS}|g;
             s|_GXFILE_|${GXFILE}|g;
             s|_GENINF_|${GENINF}|g;
             s|_PREFIX_|${OUTPRF}|g;
             s|_ST_END_|${INCSNP}|g;
             s|_TEJAAS_|${TEJAAS}|g;
             s|_SNPCUT_|${SNPCUT}|g;
             s|_GENCUT_|${GENCUT}|g;" $1 > ${JOBSUBDIR}/${JOBNAME}.bsub

        # Submit the job
        cd ${JOBSUBDIR}
        bsub < ${JOBNAME}.bsub
        cd ${CWD}
    done

done < chroms.txt
