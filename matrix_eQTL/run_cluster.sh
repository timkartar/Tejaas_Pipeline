#!/bin/sh

for CHROM in `cat chroms.txt`;
do

RUNTIME="48:00"

INPUTDIR="/cbscratch/franco/datasets/gtex"
OUTDIRBASE="/cbscratch/franco/tejaas_output/matrixEQTL_out_filteredGT_gencode/gtex_results"

# # GTFILE="${INPUTDIR}/GTEx_450Indiv_genot_imput_info04_maf01_HWEp1E6_dbSNP135IDs_donorIDs_dosage_chr${CHROM}.gz"
# GTFILE="${INPUTDIR}/prefiltered/GTEx_450Indiv_filtered_chr${CHROM}.gz"
# GXFILE="${INPUTDIR}/Whole_Blood_Analysis.v6p.normalized.expression.txt.gencode_filter"
# DONORS="${INPUTDIR}/donor_ids.fam"
# GENINF="${INPUTDIR}/genepos.gencode.v19.txt"
# DATASET="gtex"
# JOBNAME="gtex_mEQTL_${CHROM}"

# echo "Running MatrixEQTL for GTEx - CHR${CHROM}"
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



# INPUTDIR="/cbscratch/franco/datasets/cardiogenics"
# OUTDIRBASE="/cbscratch/franco/tejaas_output/matrixEQTL_out_filteredGT_gencode/cardio_results_mono_us"
# GTFILE="${INPUTDIR}/genotypes/prefiltered/CG_dosages_filtered_${CHROM}.imputed.gz"
# GXFILE="${INPUTDIR}/cardio_mono_expr.txt.gencode_filter"
# DONORS="${INPUTDIR}/genotypes/CG.sample"
# GENINF="${INPUTDIR}/../gtex/genepos.gencode.v19.txt"
# DATASET="cardiogenics"  # filtered genotypes have gtex format
# JOBNAME="cardio_mEQTL_${CHROM}"
# UDONOR="${HOME}/Tejaas_Pipeline/devtools/mono_usersamples.txt"


# echo "Running MatrixEQTL for Cardiogenics Monocytes - CHR${CHROM} "
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
#                          -i ${GENINF} \
#                          -u ${UDONOR}


INPUTDIR="/cbscratch/franco/datasets/cardiogenics"
OUTDIRBASE="/cbscratch/franco/tejaas_output/matrixEQTL_out_filteredGT_gencode/cardio_results_macro_us"
GTFILE="${INPUTDIR}/genotypes/prefiltered/CG_dosages_filtered_${CHROM}.imputed.gz"
GXFILE="${INPUTDIR}/cardio_macro_expr.txt.gencode_filter"
DONORS="${INPUTDIR}/genotypes/CG.sample"
GENINF="${INPUTDIR}/../gtex/genepos.gencode.v19.txt"
DATASET="cardiogenics"  # filtered genotypes have gtex format
JOBNAME="cardio_mEQTL_${CHROM}"
UDONOR="${HOME}/Tejaas_Pipeline/devtools/macro_usersamples.txt"

echo "Running MatrixEQTL for Cardiogenics Macrocytes- CHR${CHROM} "
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
                         -i ${GENINF} \
                         -u ${UDONOR}


done;