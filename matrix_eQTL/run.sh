
CHROM=22

INPUTDIR="${HOME}/datasets/gtex"
OUTDIRBASE="${HOME}/Tejaas_Pipeline/output/gtex_results"
GTFILE="${INPUTDIR}/GTEx_450Indiv_genot_imput_info04_maf01_HWEp1E6_dbSNP135IDs_donorIDs_dosage_chr${CHROM}.gz"
GXFILE="${INPUTDIR}/Whole_Blood_Analysis.v6p.normalized.expression.txt"
DONORS="${INPUTDIR}/donor_ids.fam"
GENINF="${INPUTDIR}/genepos.gencode.v19.txt"
DATASET="gtex"


echo "Running MatrixEQTL for GTEx "
Rscript main.R -g ${GTFILE} \
               -d ${DONORS} \
               -s ${DATASET} \
               -x ${GXFILE} \
               -o ${OUTDIRBASE} \
               -i ${GENINF}



INPUTDIR="${HOME}/datasets/cardiogenics"
OUTDIRBASE="${HOME}/Tejaas_Pipeline/output/cardio_results"
GTFILE="${INPUTDIR}/genotypes/CG_${CHROM}.imputed.gz"
GXFILE="${INPUTDIR}/cardio_mono_expr.txt"
DONORS="${INPUTDIR}/genotypes/CG.sample"
GENINF="${INPUTDIR}/../gtex/genepos.gencode.v19.txt"
DATASET="cardiogenics"

echo "Running MatrixEQTL for Cardiogenics "
Rscript main.R -g ${GTFILE} \
               -d ${DONORS} \
               -s ${DATASET} \
               -x ${GXFILE} \
               -o ${OUTDIRBASE} \
               -i ${GENINF}