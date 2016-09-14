#!/bin/bash
# Guy Leonard MMXVI
# Number of processor cores
THREADS=8

# Dependency Checks
#command -v pigz >/dev/null 2>&1 || { echo "I require pigz but it's not installed.  Aborting." >&2; exit 1;}
#command -v blastn >/dev/null 2>&1 || { echo "I require BLASTn but it's not installed.  Aborting." >&2; exit 1;}
#command -v multiqc >/dev/null 2>&1 || { echo "I require MultiQC but it's not installed.  Aborting." >&2; exit 1;}



# Working Directory
WD=$1
echo "Working Directory: $WD"

# Get filenames for current Single Cell Library
# Locations of FASTQs = Sample_**_***/raw_illumina_reads/
for DIRS in $WD/*; do

  if [ -d ${DIRS} ]; then
    echo "Working in ${DIRS}"

  ## mkdir $DIRS/raw_illumina_reads/GENES

  ## CEGMA
  # Already run, files are in
  # $DIRS/raw_illumina_reads/CEGMA

  ## SNAP 1
  # $DIRS/raw_illumina_reads/GENES/SNAP1

  ## GeneMark
  # $DIRS/raw_illumina_reads/GENES/GENEMARK

  ## MAKER 1
  # $DIRS/raw_illumina_reads/GENES/MAKER1

  ## SNAP 2
  # $DIRS/raw_illumina_reads/GENES/SNAP2

  ## AUGUSTUS
  # $DIRS/raw_illumina_reads/GENES/AUGUSTUS

  ## MAKER 2
  # $DIRS/raw_illumina_reads/GENES/MAKER2

  ## Collate GFF3 + FASTA

  fi
done
