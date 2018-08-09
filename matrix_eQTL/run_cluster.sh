#!/bin/sh

for CHROM in `cat chroms.txt`;
do

RUNTIME="48:00"

INPUTDIR="${HOME}/datasets/gtex"
# OUTDIRBASE="${HOME}/Tejaas_Pipeline/output/gtex_results"
OUTDIRBASE="/cbscratch/franco/tejaas_output/matrixEQTL_out/gtex_results"

GTFILE="${INPUTDIR}/GTEx_450Indiv_genot_imput_info04_maf01_HWEp1E6_dbSNP135IDs_donorIDs_dosage_chr${CHROM}.gz"
GXFILE="${INPUTDIR}/Whole_Blood_Analysis.v6p.normalized.expression.txt"
DONORS="${INPUTDIR}/donor_ids.fam"
GENINF="${INPUTDIR}/genepos.gencode.v19.txt"
DATASET="gtex"
JOBNAME="gtex_mEQTL_${CHROM}"

echo "Running MatrixEQTL for GTEx - CHR${CHROM}"
bsub -n 8 -a openmp -q mpi -W ${RUNTIME} -R span[hosts=1] -R cbscratch\
          -J ${JOBNAME} \
          -o ${OUTDIRBASE}/${JOBNAME}.log \
          -e ${OUTDIRBASE}/${JOBNAME}.err \
          Rscript main.R -g ${GTFILE} \
                         -c ${CHROM}  \
                         -d ${DONORS} \
                         -s ${DATASET} \
                         -x ${GXFILE} \
                         -o ${OUTDIRBASE} \
                         -i ${GENINF}



# INPUTDIR="${HOME}/datasets/cardiogenics"
# # OUTDIRBASE="${HOME}/Tejaas_Pipeline/output/cardio_results"
# OUTDIRBASE="/cbscratch/franco/tejaas_output/matrixEQTL_out/cardio_results"
# GTFILE="${INPUTDIR}/genotypes/CG_${CHROM}.imputed.gz"
# GXFILE="${INPUTDIR}/cardio_mono_expr.txt"
# DONORS="${INPUTDIR}/genotypes/CG.sample"
# GENINF="${INPUTDIR}/../gtex/genepos.gencode.v19.txt"
# DATASET="cardiogenics"
# JOBNAME="cardio_mEQTL_${CHROM}"


# echo "Running MatrixEQTL for Cardiogenics - CHR${CHROM} "
# bsub -n 8 -a openmp -q mpi -W ${RUNTIME} -R span[hosts=1] -R cbscratch\
#           -J ${JOBNAME} \
#           -o ${OUTDIRBASE}/${JOBNAME}.log \
#           -e ${OUTDIRBASE}/${JOBNAME}.err \
#           Rscript main.R -g ${GTFILE} \
#                          -c ${CHROM}  \
#                          -d ${DONORS} \
#                          -s ${DATASET} \
#                          -x ${GXFILE} \
#                          -o ${OUTDIRBASE} \
#                          -i ${GENINF}

done;