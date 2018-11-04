#!/bin/bash

source "./CONFIG.sh"

mkdir -p $COVDIR
mkdir -p $PEEREXPRDIR

while IFS='' read -r line || [ -n "$line" ]; do
    fullname=$(echo "$line" | cut -f 1 )
    shortname=$(echo "$line" | cut -f 2 )
    base=$(echo $fullname | sed 's/ - /_/g' | sed 's/ /_/g' | sed 's/(//g' | sed 's/)//g' )
    GTEXCOV="${GTEXCOVDIR}/${base}_Analysis.v6p.covariates.txt"
    COVARS="${GTEXCOVDIR}/${shortname}_nopeer_covariates.txt"

    if [ -e $GTEXCOV ]; then
        echo "Processing Covariates for Tissue: $fullname"
        # Selects only genotype PCs, sex and platform
        grep -v -i "inferred" $GTEXCOV > $COVARS

        if [ $RUNPEER = true ]; then
            # get PEER covariates correcting for covariates above
            PEERPREFIX="${shortname}_gtex.${NCOV}peer.covariates.txt"
            EXPRFILE="${EXPROUTDIR}/gtex.normalized.expression.${shortname}.txt"
            
            echo Rscript PEER.R $EXPRFILE $PEERPREFIX --n $NCOV --covar $COVARS -o $COVDIR
        fi

        # PEERresiduals="${COVDIR}/${PEERPREFIX}_PEER_residuals.txt"
        # PEEREXPR="${PEEREXPRDIR}/${shortname}_gtex.normalized.${NCOV}peer.expression.txt"
        # mv $PEERresiduals $PEEREXPR
    fi
done < $TISSUEFILE