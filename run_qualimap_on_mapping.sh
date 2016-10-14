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

    	MAPPING_DIR="$DIRS/$SAMPLE_NAME/raw_illumina_reads/BLOBTOOLS/MAPPING"
    	echo -e "\t\tBD: ${BLAST_DIR}"

	qualimap bamqc -bam scaffolds_mapped_all_reads.bam -outfile ${MAPPING_DIR}/${SAMPLE_NAME}\_qualimap.pdf
      fi
    done
  fi
done
