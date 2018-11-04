#!/bin/sh

for CHROM in `cat chroms.txt`;
do

RUNTIME="48:00"

run_cardio=false
run_us=false
run_gtex=true

wb=false
hlv=true
ms=true


if [ $run_gtex = true ]; then
     if [ $wb = true ]; then
          DATASET="gtex"
          INPUTDIR="/cbscratch/franco/datasets/gtex"
          OUTDIRPATH="/cbscratch/franco/tejaas_output/matrixEQTL_out_filteredGT_gencode/"
          OUTDIRBASE="${OUTDIRPATH}/${DATASET}_wb"

          # GTFILE="${INPUTDIR}/GTEx_450Indiv_genot_imput_info04_maf01_HWEp1E6_dbSNP135IDs_donorIDs_dosage_chr${CHROM}.gz"
          GTFILE="${INPUTDIR}/prefiltered/GTEx_450Indiv_filtered_chr${CHROM}.gz"
          GXFILE="${INPUTDIR}/Whole_Blood_Analysis.v6p.normalized.expression.txt.gencode_filter"
          DONORS="${INPUTDIR}/donor_ids.fam"
          GENINF="${INPUTDIR}/genepos.gencode.v19.txt"
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
     fi
     if [ $hlv = true ]; then
          echo "GTEx - Heart Left Ventricle"
          DATASET="gtex"
          INPUTDIR="/cbscratch/franco/datasets/gtex"
          OUTDIRPATH="/cbscratch/franco/tejaas_output/matrixEQTL_out_filteredGT_gencode/"
          OUTDIRBASE="${OUTDIRPATH}/${DATASET}_hlv"

          # GTFILE="${INPUTDIR}/GTEx_450Indiv_genot_imput_info04_maf01_HWEp1E6_dbSNP135IDs_donorIDs_dosage_chr${CHROM}.gz"
          GTFILE="${INPUTDIR}/prefiltered/GTEx_450Indiv_filtered_chr${CHROM}.gz"
          GXFILE="${INPUTDIR}/expression/GTEx_Analysis_v6p_eQTL_expression_matrices/Heart_Left_Ventricle_Analysis.v6p.normalized.expression.bed.tejaas.txt.gencode_filter"
          DONORS="${INPUTDIR}/donor_ids.fam"
          GENINF="${INPUTDIR}/genepos.gencode.v19.txt"
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
     fi

     if [ $ms = true ]; then
          echo "GTEx - Muscle Skeletal"
          DATASET="gtex"
          INPUTDIR="/cbscratch/franco/datasets/gtex"
          OUTDIRPATH="/cbscratch/franco/tejaas_output/matrixEQTL_out_filteredGT_gencode/"
          OUTDIRBASE="${OUTDIRPATH}/${DATASET}_ms"

          # GTFILE="${INPUTDIR}/GTEx_450Indiv_genot_imput_info04_maf01_HWEp1E6_dbSNP135IDs_donorIDs_dosage_chr${CHROM}.gz"
          GTFILE="${INPUTDIR}/prefiltered/GTEx_450Indiv_filtered_chr${CHROM}.gz"
          GXFILE="${INPUTDIR}/expression/GTEx_Analysis_v6p_eQTL_expression_matrices/Muscle_Skeletal_Analysis.v6p.normalized.expression.bed.tejaas.txt.gencode_filter"
          DONORS="${INPUTDIR}/donor_ids.fam"
          GENINF="${INPUTDIR}/genepos.gencode.v19.txt"
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
     fi

fi



###################################################
############# MONOCYTES ###########################
###################################################

if [ $run_cardio = true ]; then
     if [ $mono = true ]; then
          DATASET="cardiogenics"  # filtered genotypes have gtex format
          INPUTDIR="/cbscratch/franco/datasets/cardiogenics"
          OUTDIRPATH="/cbscratch/franco/tejaas_output/matrixEQTL_out_filteredGT_gencode/"
          OUTDIRBASE="${OUTDIRPATH}/${DATASET}_mono"
          GTFILE="${INPUTDIR}/genotypes/prefiltered_dosages/CG_dosages_filtered_${CHROM}.imputed.gz"
          # GXFILE="${INPUTDIR}/cardio_mono_expr.txt.gencode_filter"
          GXFILE="${INPUTDIR}/cardio_mono_matrix_nodup.txt.gencode_filter"
          DONORS="${INPUTDIR}/genotypes/CG.sample"
          GENINF="${INPUTDIR}/../gtex/genepos.gencode.v19.txt"
          JOBNAME="cardio_mEQTL_${CHROM}"


          echo "Running MatrixEQTL for Cardiogenics Monocytes - CHR${CHROM} "
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

          if [ $run_us = true ]; then
               # ~~ with user_samples ~~

               DATASET="cardiogenics"  # filtered genotypes have gtex format
               INPUTDIR="/cbscratch/franco/datasets/cardiogenics"
               OUTDIRPATH="/cbscratch/franco/tejaas_output/matrixEQTL_out_filteredGT_gencode/"
               OUTDIRBASE="${OUTDIRPATH}/${DATASET}_mono_us"
               GTFILE="${INPUTDIR}/genotypes/prefiltered_dosages/CG_dosages_filtered_${CHROM}.imputed.gz"
               # GXFILE="${INPUTDIR}/cardio_mono_expr.txt.gencode_filter"
               GXFILE="${INPUTDIR}/cardio_mono_matrix_nodup.txt.gencode_filter"
               DONORS="${INPUTDIR}/genotypes/CG.sample"
               GENINF="${INPUTDIR}/../gtex/genepos.gencode.v19.txt"
               JOBNAME="cardio_mEQTL_${CHROM}"
               UDONOR="${HOME}/Tejaas_Pipeline/devtools/mono_usersamples.txt"


               echo "Running MatrixEQTL for Cardiogenics Monocytes - CHR${CHROM} "
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
          fi
     fi

     if [ $macro = true ]; then
          ###################################################
          ############# MACROPHAGES #########################
          ###################################################

          DATASET="cardiogenics"  # filtered genotypes have gtex format
          INPUTDIR="/cbscratch/franco/datasets/cardiogenics"
          OUTDIRPATH="/cbscratch/franco/tejaas_output/matrixEQTL_out_filteredGT_gencode/"
          OUTDIRBASE="${OUTDIRPATH}/${DATASET}_macro"
          GTFILE="${INPUTDIR}/genotypes/prefiltered_dosages/CG_dosages_filtered_${CHROM}.imputed.gz"
          # GXFILE="${INPUTDIR}/cardio_macro_expr.txt.gencode_filter"
          GXFILE="${INPUTDIR}/cardio_macro_matrix_nodup.txt.gencode_filter"
          DONORS="${INPUTDIR}/genotypes/CG.sample"
          GENINF="${INPUTDIR}/../gtex/genepos.gencode.v19.txt"
          JOBNAME="cardio_mEQTL_${CHROM}"

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
                                   -i ${GENINF} 


          # ~~ with user_samples ~~
          if [ $run_us = true ]; then
               DATASET="cardiogenics"  # filtered genotypes have gtex format
               INPUTDIR="/cbscratch/franco/datasets/cardiogenics"
               OUTDIRPATH="/cbscratch/franco/tejaas_output/matrixEQTL_out_filteredGT_gencode/"
               OUTDIRBASE="${OUTDIRPATH}/${DATASET}_macro_us"
               GTFILE="${INPUTDIR}/genotypes/prefiltered_dosages/CG_dosages_filtered_${CHROM}.imputed.gz"
               # GXFILE="${INPUTDIR}/cardio_macro_expr.txt.gencode_filter"
               GXFILE="${INPUTDIR}/cardio_macro_matrix_nodup.txt.gencode_filter"
               DONORS="${INPUTDIR}/genotypes/CG.sample"
               GENINF="${INPUTDIR}/../gtex/genepos.gencode.v19.txt"
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
          fi
     fi
fi

done;