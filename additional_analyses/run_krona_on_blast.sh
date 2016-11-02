#!/bin/bash
# Guy Leonard MMXVI
# Number of processor cores

# Working Directory
WD=`pwd`
echo "Working Directory: $WD"
# Script Dir
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Get dirnames for current Single Cell Library
# Locations of FASTQs = Sample_**_***/raw_illumina_reads/
for DIRS in $WD/*; do
  if [ -d ${DIRS} ]; then
    echo "Working in ${DIRS}"
    GROUP_NAME="$(basename $DIRS)"
    echo "GN: $GROUP_NAME"

    for SAMPLES in $GROUP_NAME/*; do
      if [ -d ${SAMPLES} ]; then 
        SAMPLE_NAME="$(basename $SAMPLES)"
  	echo -e "\tSample: $SAMPLE_NAME"

    	BLAST_DIR="$DIRS/$SAMPLE_NAME/raw_illumina_reads/BLOBTOOLS/BLAST"
    	echo -e "\t\tBD: ${BLAST_DIR}"

        cat ${BLAST_DIR}/scaffolds_vs_nt_1e-10.megablast | cut -d$'\t' -f1,2,13 > ${BLAST_DIR}/${SAMPLE_NAME}\_krona.tsv

        ktImportTaxonomy ${BLAST_DIR}/${SAMPLE_NAME}\_krona.tsv -i -o ${BLAST_DIR}/${SAMPLE_NAME}.html
      fi
    done
  fi
done
