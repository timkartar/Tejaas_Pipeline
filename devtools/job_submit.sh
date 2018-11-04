#!/bin/bash

beta=$2
METHOD=$3
INPUTDIR="${INPUTDIR}"
OUTDIRBASE="${OUTDIRBASE}/${METHOD}/beta_${beta}"
JOBSUBDIR="${OUTDIRBASE}/jobsub"
CWD=`pwd`

MAX_NSNP=100000
SNPCUT=0.001
GENCUT=0.001

if [ ! -d ${JOBSUBDIR} ]; then
    mkdir -p ${JOBSUBDIR}
fi
    
while read j; do

    echo "Submitting jobs for Chromosome ${j}."

    CHROM=$j
    GTFILE="${GTPATH/__CHROM__/$CHROM}"
    GXFILE="${GXFILE}"
    DONORS="${DONORS}"
    GENINF="${GENINF}"
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
             s|_GENCUT_|${GENCUT}|g;
             s|_METHOD_|${METHOD}|g;
             s|_CHROM_|${CHROM}|g;
             s|_XTRAF_|${EXTRAFLAGS}|g;
             s|_OPTSAMPL_|${OPT_SAMPLES}|g;
             s|_BETA_|${beta}|g;" $1 > ${JOBSUBDIR}/${JOBNAME}.bsub

        # Submit the job
        cd ${JOBSUBDIR}
        ok=0
        # workaround, sometimes the file was submitted to the cluster but not written to disk yet
        while [ $ok -lt 1 ]; do
            if [ -s ${JOBNAME}.bsub ]; then
                bsub < ${JOBNAME}.bsub
                ok=1
            else
                echo "File is empty, retrying.."
                sleep 1;
            fi
        done
        cd ${CWD}
    done

done < chroms.txt
