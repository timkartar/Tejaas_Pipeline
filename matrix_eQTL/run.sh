
CHROM=6

# for CHROM in `cat ../devtools/chroms.txt`;
# do

# INPUTDIR="/cbscratch/franco/datasets/gtex"
# OUTDIRBASE="/cbscratch/franco/tejaas_output/matrixEQTL_out/gtex_results"
# GTFILE="${INPUTDIR}/GTEx_450Indiv_genot_imput_info04_maf01_HWEp1E6_dbSNP135IDs_donorIDs_dosage_chr${CHROM}.gz"
# GXFILE="${INPUTDIR}/Whole_Blood_Analysis.v6p.normalized.expression.txt"
# DONORS="${INPUTDIR}/donor_ids.fam"
# GENINF="${INPUTDIR}/genepos.gencode.v19.txt"
# DATASET="gtex"


# echo "Running MatrixEQTL for GTEx - CHR${CHROM}"
# Rscript main.R -g ${GTFILE} \
#                -c ${CHROM}  \
#                -d ${DONORS} \
#                -s ${DATASET} \
#                -x ${GXFILE} \
#                -o ${OUTDIRBASE} \
#                -i ${GENINF}



INPUTDIR="/cbscratch/franco/datasets/cardiogenics"

# OUTDIRBASE="/cbscratch/franco/tejaas_output/matrixEQTL_out/cardio_results"
OUTDIRBASE="/cbscratch/franco/tejaas_output/matrixEQTL_out/cardio_results_macro_us"

GTFILE="${INPUTDIR}/genotypes/prefiltered_dosages/CG_${CHROM}.imputed.gz"

# GXFILE="${INPUTDIR}/cardio_mono_expr.txt"
GXFILE="${INPUTDIR}/cardio_macro_expr.txt"

DONORS="${INPUTDIR}/genotypes/CG.sample"
GENINF="${INPUTDIR}/../gtex/genepos.gencode.v19.txt"
DATASET="cardiogenics"
UDONOR="${HOME}/Tejaas_Pipeline/devtools/usersamples_macro.txt"



echo "Running MatrixEQTL for Cardiogenics - CHR${CHROM} "
Rscript main.R -g ${GTFILE} \
               -c ${CHROM}  \
               -d ${DONORS} \
               -s ${DATASET} \
               -x ${GXFILE} \
               -o ${OUTDIRBASE} \
               -i ${GENINF} \
               -u ${UDONOR}

# done;