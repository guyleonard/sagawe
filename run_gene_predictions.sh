#!/bin/bash
# Guy Leonard MMXVI
# Number of processor cores
THREADS=8

# Dependency Checks
command -v pigz >/dev/null 2>&1 || { echo "I require pigz but it's not installed.  Aborting." >&2; exit 1;}

command -v blastn >/dev/null 2>&1 || { echo "I require BLASTn but it's not installed.  Aborting." >&2; exit 1;}

command -v multiqc >/dev/null 2>&1 || { echo "I require MultiQC but it's not installed.  Aborting." >&2; exit 1;}



# Working Directory
WD=$1
echo "$WD"

# Get filenames for current Single Cell Library
# Locations of FASTQs = Sample_**_***/raw_illumina_reads/
for DIRS in $WD ; do
        echo "Working in $DIRS"
done
