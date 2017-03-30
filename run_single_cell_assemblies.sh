#!/usr/bin/env bash

############################
## User Defined Variables ##
############################

# NCBI 'nt' Database Location and name (no extension)
NCBI_NT="/home/ubuntu/blast/nt"
# NCBI Taxonomy Dump
# ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
NCBI_TAXDMP="/home/ubuntu/blast/taxonomy"
# NCBI Taxonomy DB Location
# ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
NCBI_TAXDB="/home/ubuntu/blast/taxdb"
export BLASTDB=$NCBI_TAXDB
# CEGMA DIR
CEGMA_DIR="/home/ubuntu/single_cell_workflow/build/CEGMA_v2.5"
# BUSCO Lineage Location
BUSCO_V1_DB="/home/ubuntu/busco/v1"
BUSCO_V2_DB="/home/ubuntu/busco/v2"
# Augustus Config Path
AUGUSTUS_CONFIG_PATH="/home/ubuntu/single_cell_workflow/build/augustus-3.2.2/config"


#############################################
## Program Functions - Do Not Change Below ##
##         Roughly In Order of Use         ##
#############################################

## Run bbnorm
# This is to digitally normalise your
# data if it will not assemble, remember
# you should map the original reads to an
# assembly and not the normalised reads,
# this script does not currently do that
# so as to to preserve overlap/trimmed reads
# used in the assembly for normal analysis...
function run_normalisation () {
    normalised_dir="$output_dir/normalised"

    if [ -d $normalised_dir ] ; then
        echo "Normalisation Previously Run, next..."
        NORMALISED="true"
    else
        mkdir -p "$normalised_dir"

        echo "Running BBNorm"
        bbnorm.sh in1="$READ1" in2="$READ2" \
        out1="$normalised_dir/${READ1/.gz/.norm.gz}" \
        out2="$normalised_dir/${READ2/.gz/.norm.gz}" \
        outt="$normalised_dir/excluded_reads.fastq.gz" \
        hist="$normalised_dir/input_kmer_depth.hist" \
        histout="$normalised_dir/output_kmer_depth.hist" \
        threads="$THREADS"
    fi
}

## Function to control which overlapper is used
function run_overlapper () {
    if [ $overlap_option == "bbmerge" ] ; then
        run_bbmerge
    elif [ $overlap_option == "pear" ] ; then
        run_pear
    else
        run_pear
    fi
}

## Run Pear Overlapper
# default settings
function run_pear () {
    overlapped_dir="$output_dir/overlapped"
    mkdir -p "$overlapped_dir"

    absolute_path="$( cd "$overlapped_dir" && pwd )"

    if [ "$normalised" == 'true' ] ; then
        echo "Running PEAR Assembler with Normalised Reads"
        pear -f "$normalised_dir/${READ1/.gz/.norm.gz}" -r "$normalised_dir/${READ2/.gz/.norm.gz}" \
        -o "$overlapped_dir/pear_overlap" -j "$THREADS" | tee "$overlapped_dir/pear.log"
    else
        echo "Running PEAR Assembler"
        pear -f "$READ1" -r "$READ2" \
        -o "$overlapped_dir/pear_overlap" -j "$THREADS" | tee "$overlapped_dir/pear.log"
    fi

    gzip_fastq

    ln -s "$absolute_path/pear_overlap.assembled.fastq.gz" "$absolute_path/assembled.fastq.gz"
    ln -s "$absolute_path/pear_overlap.unassembled.forward.fastq.gz" "$absolute_path/unassembled.forward.fastq.gz"
    ln -s "$absolute_path/pear_overlap.unassembled.reverse.fastq.gz" "$absolute_path/unassembled.reverse.fastq.gz"
}

## Run BBMerge
# Roughly same settings as PEAR
function run_bbmerge () {
    overlapped_dir="$output_dir/overlapped"
    mkdir -p "$overlapped_dir"

    absolute_path="$( cd "$overlapped_dir" && pwd )"

    if [ "$normalised" == 'true' ] ; then
        echo "Running BBMerge with Normalised Reads"
        bbmerge.sh in1="$normalised_dir/${READ1/.gz/.norm.gz}" in2="$normalised_dir/${READ2/.gz/.norm.gz}" \
        out="$absolute_path/assembled.fastq.gz" \
        outu1="$absolute_path/unassembled.forward.fastq.gz" \
        outu2="$absolute_path/unassembled.reverse.fastq.gz" \
        minoverlap=10 ziplevel=9
    else
        echo "Running BBMerge"
        bbmerge.sh in1="$READ1" in2="$READ2" \
        out="$absolute_path/assembled.fastq.gz" \
        outu1="$absolute_path/unassembled.forward.fastq.gz" \
        outu2="$absolute_path/unassembled.reverse.fastq.gz" \
        minoverlap=10 ziplevel=9
    fi
}

function gzip_fastq () {
    overlapped_dir="$output_dir/overlapped"
    absolute_path="$( cd "$overlapped_dir" && pwd )"

    if command -v "clumpify.sh" >/dev/null 2>&1 ; then
        echo "Running clumpify.sh"
        clumpify.sh in="$absolute_path/pear_overlap.assembled.fastq" \
        out="$absolute_path/pear_overlap.assembled.fastq.gz"

        clumpify.sh in="$absolute_path/pear_overlap.unassembled.forward.fastq" \
        in2="$absolute_path/pear_overlap.unassembled.reverse.fastq" \
        out="$absolute_path/pear_overlap.unassembled.forward.fastq.gz" \
        out2="$absolute_path/pear_overlap.unassembled.reverse.fastq.gz"
    elif command -v "pigz" >/dev/null 2>&1 ; then
        echo "Running pigz"
        pigz -9 -R "$absolute_path/pear_overlap.assembled.fastq"
        pigz -9 -R "$absolute_path/pear_overlap.unassembled.forward.fastq"
        pigz -9 -R "$absolute_path/pear_overlap.unassembled.reverse.fastq"
    else
        echo "Running gzip"
        gzip -9 --rsyncable "$absolute_path/pear_overlap.assembled.fastq"
        gzip -9 --rsyncable "$absolute_path/pear_overlap.unassembled.forward.fastq"
        gzip -9 --rsyncable "$absolute_path/pear_overlap.unassembled.reverse.fastq"
    fi
}

## Run Trim Galore!
# Two sets of trimming:
# Assembled/overlapped reads
# Unassembled paired reads.
# Minimum length of 100
# Minimum quality of Q20
# run FASTQC on trimmed
# GZIP output
function run_trim_galore () {
    overlapped_dir="$output_dir/overlapped"

    if [ ! -f "$overlapped_dir/assembled.fastq.gz" ] ; then
        echo "Read Overlapper Was/Did Not Run. Please use -p option."
    else
        trimmed_dir="$output_dir/trimmed"
        mkdir -p "$trimmed_dir"
        fastqc_dir="$trimmed_dir/fastqc"
        mkdir -p "$fastqc_dir"

        echo "Running Trim Galore! on Untrimmed Assembled Reads"
        trim_galore -q 20 --fastqc --'gzip' --length 100 \
        "$overlapped_dir/assembled.fastq.gz" -o "$trimmed_dir" | tee "$trimmed_dir/trim_galore_assembled.log"

        echo "Running Trim Galore! on Untrimmed Un-assembled Reads"
        trim_galore -q 20 --fastqc --'gzip' --length 100 --paired --retain_unpaired \
        "$overlapped_dir/unassembled.forward.fastq.gz" "$overlapped_dir/unassembled.reverse.fastq.gz" \
        -o "$trimmed_dir" | tee "$trimmed_dir/trim_galore_unassembled.log"

        absolute_path="$( cd "$trimmed_dir" && pwd )"

        mv "$absolute_path/"*".html" "$absolute_path/fastqc"
        mv "$absolute_path/"*".zip" "$absolute_path/fastqc"
        mv "$absolute_path/"*".txt" "$absolute_path/fastqc"

        # for spades, we should include the now unpaired reads as one library
        # it makes sense but at the same time freezes my brain
        # you can cat .gz files together!
        cat "$trimmed_dir/unassembled.forward_unpaired_1.fq.gz" \
        "$trimmed_dir/unassembled.reverse_unpaired_2.fq.gz" > "$trimmed_dir/unassembled_unpaired.fq.gz"
    fi
}

## Run SPAdes
# single cell mode - default kmers 21,33,55
# careful mode - runs mismatch corrector
# use all 4 sets of output reads
# 1 x overlapped
# 2 x unoverlapped paired
# 1 x unpaired
function run_assembly_spades () {
    trimmed_dir="$output_dir/trimmed"

    if [ ! -f "$trimmed_dir/assembled_trimmed.fq.gz" ] ; then
        echo "Read Overlapper Was/Did Not Run. Please use -p option."
    else
        assembly_dir="$output_dir/assembly"
        mkdir -p "$assembly_dir"
        echo "Running SPAdes"
        spades.py --sc --careful -t "$THREADS" \
        --s1 "$trimmed_dir/assembled_trimmed.fq.gz" \
        --pe1-1 "$trimmed_dir/unassembled.forward_val_1.fq.gz" \
        --pe1-2 "$trimmed_dir/unassembled.reverse_val_2.fq.gz" \
        --pe1-s "$trimmed_dir/unassembled_unpaired.fq.gz" \
        -o "$assembly_dir"
    fi
}

# Run QUAST
# eukaryote mode
# glimmer protein predictions
function report_quast () {
    assembly_dir="$output_dir/assembly"

    if [ ! -f "$assembly_dir/$ASSEMBLY" ] ; then
        echo -e "[ERROR]: SPAdes scaffolds cannot be found. Aborting."
        exit 1
    else
        quast_dir="$output_dir/reports/quast"
        mkdir -p "$quast_dir"
        echo "Running QUAST"
        quast.py -o "$quast_dir" -t "$THREADS" \
        --min-contig 100 -f --eukaryote $QUAST_SCAFFOLDS \
        --glimmer "$assembly_dir/$ASSEMBLY" | tee "$quast_dir/quast.log"
    fi
}

# Run CEGMA
# Genome mode
function report_cegma () {
    assembly_dir="$output_dir/assembly"

    if [ ! -f "$assembly_dir/$ASSEMBLY" ] ; then
        echo -e "[ERROR]: SPAdes scaffolds cannot be found. Aborting."
        exit 1
    else
        cegma_dir="$output_dir/reports/cegma"
        mkdir -p "$cegma_dir"

        echo "Running CEGMA"
        cegma -T "$THREADS" -g "$assembly_dir/$ASSEMBLY" -o "$cegma_dir/cegma" | tee "$cegma_dir/cegma.log"

        # Tidy up, CEGMA's -o option doesn't seem to output to a dir!?
        absolute_path="$( cd "$output_dir" && pwd )"
        mv "$absolute_path/cegma."* "$absolute_path/reports/cegma"
    fi
}

# Run BUSCO v1 - LEGACY
function report_busco_v1 () {
    assembly_dir="$output_dir/assembly"

    if [ ! -f "$assembly_dir/$ASSEMBLY" ] ; then
        echo -e "[ERROR]: SPAdes scaffolds cannot be found. Aborting."
        exit 1
    else
        busco_dir="$output_dir/reports/busco"
        mkdir -p "$busco_dir"

        absolute_path="$( cd "$output_dir" && pwd )"
        current_dir=$(pwd)

        IFS=\, read -a current_db <<<"$busco_v2_dbs"

        for x in "${current_db[@]}";do

            BUSCO_v1.22.py \
            -g "$assembly_dir/$ASSEMBLY" \
            -c "$THREADS" -l "$BUSCO_V1_DB/$y" -o "$y" -f | tee "$busco_dir/busco.log"

            # Tidy up busco
            mv "$current_dir/run_$y" "$absolute_path/reports/busco"
        done
    fi
}

# Run BUSCO v2
function report_busco_v2 () {
    assembly_dir="$output_dir/assembly"

    if [ ! -f "$assembly_dir/$ASSEMBLY" ] ; then
        echo -e "[ERROR]: SPAdes scaffolds cannot be found. Aborting."
        exit 1
    else
        busco_dir="$output_dir/reports/busco"
        mkdir -p "$busco_dir"
        mkdir -p "$busco_dir/summaries"

        absolute_path="$( cd "$output_dir" && pwd )"
        current_dir=$(pwd)

        # this needs a default!!! eukaryota_odb9,protists_ensembl

        IFS=\, read -a current_db <<<"$busco_v2_dbs"

        for x in "${current_db[@]}";do

            BUSCO.py \
            -i "$assembly_dir/$ASSEMBLY" -m genome \
            -c "$THREADS" -l "$BUSCO_V2_DB/$x" -o "$x" -f | tee "$busco_dir/busco_$x.log"

            # Tidy up busco, it won't take a path as an output, so
            # resorts to the basename dir of the current dir the script
            # is run from for out in a "run_" folder - somewhat annoying.
            mv "$current_dir/run_$x" "$absolute_path/reports/busco"

            ln -s "$absolute_path/reports/busco/run_$x/short_summary_$x.txt" "$busco_dir/summaries"
        done

        BUSCO_plot.py \
        -wd "$busco_dir/summaries"
    fi
}

function report_multiqc () {
    multiqc -o "$output_dir/reports" -z "$output_dir"
}

function blobtools_bwa () {
    assembly_dir="$output_dir/assembly"

    if [ ! -f "$assembly_dir/$ASSEMBLY" ] ; then
        echo -e "[ERROR]: SPAdes scaffolds cannot be found. Aborting."
        exit 1
    else
        blobtools_dir="$output_dir/reports/blobtools"
        mkdir -p "$blobtools_dir"

        blobtools_map="$blobtools_dir/mapping"
        mkdir -p "$blobtools_map"

        absolute_path="$( cd "$output_dir" && pwd )"

        trimmed_dir="$absolute_dir/trimmed"

        ln -s "$absolute_path/assembly/$ASSEMBLY" "$absolute_path/reports/blobtools/mapping/$ASSEMBLY"

        # index assembly (scaffolds.fa) with BWA
        echo "Indexing Assembly"
        bwa index -a bwtsw "$blobtools_map/$ASSEMBLY" | tee "$blobtools_map/bwa.log"

        # map reads to assembly with BWA MEM
        # we will have to do this for all 5 sets of reads and then merge
        echo "Mapping Assembled reads to Assembly"
        bwa mem -t "$THREADS" "$blobtools_map/$ASSEMBLY" \
        "$trimmed_dir/assembled_trimmed.fq.gz" \
        > "$blobtools_map/scaffolds_mapped_assembled_reads.sam" | tee -a "$blobtools_map/bwa.log"

        echo "Mapping Un-assembled & Un-Paired reads to Assembly - Forward"
        bwa mem -t "$THREADS" "$blobtools_map/$ASSEMBLY" \
        "$trimmed_dir/unassembled.forward_unpaired_1.fq.gz" \
        > "$blobtools_map/scaffolds_mapped_unassembled_unpaired_forward_reads.sam" | tee -a "$blobtools_map/bwa.log"

        echo "Mapping Un-assembled & Un-Paired reads to Assembly - Reverse"
        bwa mem -t "$THREADS" "$blobtools_map/$ASSEMBLY" \
        "$trimmed_dir/unassembled.reverse_unpaired_2.fq.gz" \
        > "$blobtools_map/scaffolds_mapped_unassembled_unpaired_reverse_reads.sam" | tee -a "$blobtools_map/bwa.log"

        echo "Mapping Un-assembled but still Paired reads to Assembly"
        bwa mem -t "$THREADS" "$blobtools_map/$ASSEMBLY" \
        "$trimmed_dir/unassembled.forward_val_1.fq.gz" \
        "$trimmed_dir/unassembled.reverse_val_2.fq.gz" \
        > "$blobtools_map/scaffolds_mapped_unassembled_paired_reads.sam" | tee -a "$blobtools_map/bwa.log"
    fi
}

function blobtools_samtools () {
    blobtools_dir="$output_dir/reports/blobtools"
    blobtools_map="$blobtools_dir/mapping"

    # sort and convert sam to bam with SAMTOOLS
    echo "Sorting Assembled SAM File and Converting to BAM"
    samtools sort -@ "$THREADS" -o "$blobtools_map/scaffolds_mapped_assembled_reads.bam" \
    "$blobtools_map/scaffolds_mapped_assembled_reads.sam" | tee -a "$blobtools_map/samtools.log"

    # sort and convert sam to bam with SAMTOOLS
    echo "Sorting Un-assembled & Un-Paired SAM Files and Converting to BAM - Forward"
    samtools sort -@ "$THREADS" -o "$blobtools_map/scaffolds_mapped_unassembled_unpaired_forward_reads.bam" \
    "$blobtools_map/scaffolds_mapped_unassembled_unpaired_forward_reads.sam" | tee -a "$blobtools_map/samtools.log"

    echo "Sorting Un-assembled & Un-Paired SAM Files and Converting to BAM - Reverse"
    samtools sort -@ "$THREADS" -o "$blobtools_map/scaffolds_mapped_unassembled_unpaired_reverse_reads.bam" \
    "$blobtools_map/scaffolds_mapped_unassembled_unpaired_reverse_reads.sam" | tee -a "$blobtools_map/samtools.log"

    # sort and convert sam to bam with SAMTOOLS
    echo "Sorting Un-assembled & Paired SAM File and Converting to BAM"
    samtools sort -@ "$THREADS" -o "$blobtools_map/scaffolds_mapped_unassembled_paired_reads.bam" \
    "$blobtools_map/scaffolds_mapped_unassembled_paired_reads.sam" | tee -a "$blobtools_map/samtools.log"

    # Merge SAM files
    echo "Merging 4 BAM files"
    samtools merge -@ "$THREADS" -f "$blobtools_map/scaffolds_mapped_all_reads.bam" \
    "$blobtools_map/scaffolds_mapped_assembled_reads.bam" \
    "$blobtools_map/scaffolds_mapped_unassembled_paired_reads.bam" \
    "$blobtools_map/scaffolds_mapped_unassembled_unpaired_forward_reads.bam" \
    "$blobtools_map/scaffolds_mapped_unassembled_unpaired_reverse_reads.bam" | tee -a "$blobtools_map/samtools.log"

    echo "Indexing Bam"
    samtools index "$blobtools_map/scaffolds_mapped_all_reads.bam" | tee -a "$blobtools_map/samtools.log"

    # delete sam files
    rm "$blobtools_map/scaffolds_mapped_assembled_reads.sam"
    rm "$blobtools_map/scaffolds_mapped_unassembled_paired_reads.sam"
    rm "$blobtools_map/scaffolds_mapped_unassembled_unpaired_forward_reads.sam"
    rm "$blobtools_map/scaffolds_mapped_unassembled_unpaired_reverse_reads.sam"
}

function blobtools_blast () {
    blobtools_dir="$output_dir/reports/blobtools"
    blobtools_map="$blobtools_dir/mapping"

    blobtools_blast="$blobtools_dir/blast"
    mkdir -p "$blobtools_blast"

    echo "Running BLAST"
    blastn -task megablast \
    -query "$blobtools_map/$ASSEMBLY" \
    -db "$NCBI_NT/nt" \
    -evalue 1e-10 \
    -num_threads "$THREADS" \
    -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
    -culling_limit 5 \
    -out "$blobtools_blast/scaffolds_vs_nt_1e-10.megablast" | tee "$blobtools_blast/blast.log"
}

function blobtools_create () {
    blobtools_dir="$output_dir/reports/blobtools"
    blobtools_map="$blobtools_dir/mapping"
    blobtools_blast="$blobtools_dir/blast"

    echo "Running BlobTools CREATE - slow"
    blobtools create -i "$blobtools_map/$ASSEMBLY" \
    --nodes "$NCBI_TAXDMP/nodes.dmp" --names "$NCBI_TAXDMP/names.dmp" \
    -t "$blobtools_blast/scaffolds_vs_nt_1e-10.megablast" \
    -b "$blobtools_map/scaffolds_mapped_all_reads.bam" \
    -o "$blobtools_dir/scaffolds_mapped_reads_nt_1e-10_megablast_blobtools" | tee -a "$blobtools_dir/blobtools.log"
}

function blobtools_table () {
    blobtools_dir="$output_dir/reports/blobtools"
    mkdir -p "$blobtools_dir/table"

    # Standard Output - Phylum
    echo "Running BlobTools View"
    blobtools view -i "$blobtools_dir/scaffolds_mapped_reads_nt_1e-10_megablast_blobtools.blobDB.json" \
    --out "$blobtools_dir/table/scaffolds_mapped_reads_nt_1e-10_megablast_blobtools_phylum_table.csv" | tee -a "$blobtools_dir/blobtools.log"

    # Other Output - Species
    blobtools view -i "$blobtools_dir/scaffolds_mapped_reads_nt_1e-10_megablast_blobtools.blobDB.json" \
    --out "$blobtools_dir/table/scaffolds_mapped_reads_nt_1e-10_megablast_blobtools_superkingdom_table.csv" \
    --rank superkingdom | tee -a "$blobtools_dir/blobtools.log"
}

# run blobtools plot - image output
function blobtools_image () {
    blobtools_dir="$output_dir/reports/blobtools"
    mkdir -p "$blobtools_dir/images/png"
    mkdir -p "$blobtools_dir/images/svg"

    # Standard Output - Phylum, 7 Taxa
    echo "Running BlobTools Plots - Standard + SVG"
    blobtools plot -i "$blobtools_dir/scaffolds_mapped_reads_nt_1e-10_megablast_blobtools.blobDB.json" \
    -o "$blobtools_dir/images/png/" | tee -a "$blobtools_dir/blobtools.log"

    blobtools plot -i "$blobtools_dir/scaffolds_mapped_reads_nt_1e-10_megablast_blobtools.blobDB.json" \
    --format svg -o "$blobtools_dir/images/svg/" | tee -a "$blobtools_dir/blobtools.log"


    # Other Output - Super Kingdom, 7 Taxa
    echo "Running BlobTools Plots - SuperKingdom + SVG"
    blobtools plot -i "$blobtools_dir/scaffolds_mapped_reads_nt_1e-10_megablast_blobtools.blobDB.json" \
    -r superkingdom -o "$blobtools_dir/images/png/" | tee -a "$blobtools_dir/blobtools.log"

    blobtools plot -i "$blobtools_dir/scaffolds_mapped_reads_nt_1e-10_megablast_blobtools.blobDB.json" \
    -r superkingdom --format svg -o "$blobtools_dir/images/svg/" | tee -a "$blobtools_dir/blobtools.log"
}

function report_blobtools () {
    blobtools_dir="$output_dir/reports/blobtools"

    blobtools_bwa
    blobtools_samtools
    blobtools_blast
    blobtools_create
    blobtools_table
    blobtools_image

    absolute_path="$( cd "$output_dir" && cd ../ && pwd )"

    # tidy up one file blobtools creates in wrong folder
    mv "$absolute_path/scaffolds_mapped_all_reads.bam.cov" "$blobtools_dir"
}

#########################
## Accessory Functions ##
#########################

function check_exe () {
    program=$1
    command -v "$program" >/dev/null 2>&1 || { echo "$program is required, but it is not in PATH.  Please install/ammend." >&2; exit 1;}
}

function cores () {
    cores=$(grep -c ^processor /proc/cpuinfo 2>/dev/null || sysctl -n hw.ncpu)
    echo $(($cores / 2))
}

function export_cegma () {
    if [ ! -d "$CEGMA_DIR" ] ; then
        echo "[ERROR]: Incorrect CEGMA Path. Is your path correct?"
        echo "$CEGMA_DIR"
        exit 1
    else
        export CEGMA="$CEGMA_DIR"
        export PATH=$PATH:"$CEGMA_DIR"/bin
        export PERL5LIB=$PERL5LIB:"$CEGMA_DIR"/lib
        export WISECONFIGDIR=/usr/share/wise/
    fi
}

function ncbi_nt () {
    if [ ! -f "$NCBI_NT/nt.nal" ] ; then
        echo "[ERROR]: Missing NCBI NT Libraries. Is your path correct?"
        echo "$NCBI_NT"
        exit 1
    fi
}

function ncbi_taxonomy () {
    if [ ! -f "$NCBI_TAXDMP/nodes.dmp" ] ; then
        echo "[ERROR]: Missing NCBI Taxonomy Libraries. Is your path correct?"
        echo "$NCBI_TAXDMP"
        exit 1
    fi
}

function busco_db () {
    if [ ! -d "$BUSCO_V1_DB" ] ; then
        echo "[ERROR]: Missing BUSCO Lineage Directory. Is your path correct?"
        echo "$BUSCO_V1_DB"
        exit 1
    fi
}

function augustus () {
    if [ ! -d "$AUGUSTUS_CONFIG_PATH" ] ; then
        echo "[ERROR]: Missing AUGUSTUS_CONFIG_PATH Directory. Is your path correct?"
        echo "$AUGUSTUS_CONFIG_PATH"
        exit 1
    else
        export AUGUSTUS_CONFIG_PATH="$AUGUSTUS_CONFIG_PATH"
    fi
}

function help_message () {
    echo -e "Single Amplified Genome Assembly Workflow"
    echo -e "Required Option:"
    echo -e "  -o <output_dir>\tOutput Directory"
    echo -e "File Options:"
    echo -e "  -f <forward.fastq>\tForward Reads"
    echo -e "  -r <reverse.fastq>\tReverse Reads"
    echo -e "Optional Parameters (ordered):"
    echo -e "  -n \tRead Normalisation"
    echo -e "  -S \tUse Scaffolds Instead of Contigs"
    echo -e "  -a \tRun All Options Below (p{pear}tsqcb{eukaryota_odb9}Bm)"
    echo -e "Workflow Parameters:"
    echo -e "  -p <pear|bbmerge>\tOverlap Reads"
    echo -e "  -t \tTrim Overlapped Reads"
    echo -e "  -s \tAssemble Trimmed Reads"
    echo -e "Reports:"
    echo -e "  -q \tRun QUAST"
    echo -e "  -c \tRun CEGMA"
    echo -e "  -b <db1,db2,...>\tRun BUSCO v2"
    echo -e "  -B \tRun BlobTools"
    echo -e "  -m \tRun MultiQC"
    echo -e "Legacy (soon to be deprecated):"
    echo -e "  -l \tRun BUSCO v1 - legacy"
    echo -e "Example: run_single_cell_assemblies.sh -f r1.fastq -r r2.fastq -o output_dir -n -S -a"
    exit 1
}

##############################
## Initial Workflow Actions ##
##############################

# Set Cores
THREADS=$(cores)

# Check that we have the required programs
exes=(pigz clumpify.sh trim_galore pear spades.py quast.py cegma BUSCO_v1.22.py BUSCO.py BUSCO_plot.py bwa samtools blastn blobtools multiqc)
for program in ${exes[@]} ; do
    check_exe "$program"
done

NORMALISED="false"

# Check for correct Paths and Exports
export_cegma
ncbi_taxonomy
ncbi_nt
busco_db
augustus

NUMARGS=$#
if [ "$NUMARGS" -eq 0 ]; then
  help_message
fi

# defaulting to using contigs over scaffolds
# mainly to remove the Ns as inflation from 
# statistics
ASSEMBLY="contigs.fasta"

######################
## Pipeline Options ##
######################
while getopts f:r:o:np:tsSqclb:Bmah FLAG; do
    case $FLAG in
        f)
            READ1=$OPTARG
            ;;
        r)
            READ2=$OPTARG
            ;;
        o)
            output_dir=$OPTARG
            mkdir -p "$output_dir"
            echo "Single Amplified Genome Assembly Workflow" > $output_dir/$output_dir.log
            echo -e "\thttps://github.com/guyleonard/single_cell_workflow\n" >> $output_dir/$output_dir.log
            echo "Run: $(date +%F+%R)" >> $output_dir/$output_dir.log
            echo "Command: ${0} ${@}" >> $output_dir/$output_dir.log
            ;;
        n)
            run_normalisation
            NORMALISED="true"
            ;;
        p)
            overlap_option=$OPTARG
            run_overlapper
            ;;
        t)
            run_trim_galore
            ;;
        s)
            run_assembly_spades
            ;;
        q)
            report_quast
            ;;
        c)
            report_cegma
            ;;
        l)
            busco_v1_dbs=$OPTARG
            report_busco_v1
            ;;
        b)
            busco_v2_dbs=$OPTARG
            report_busco_v2
            ;;
        B)
            report_blobtools
            ;;
        m)
            report_multiqc
            ;;
        S)
            ASSEMBLY="scaffolds.fasta"
            QUAST_SCAFFOLDS="--scaffolds"
            ;;
        a)
             overlap_option="pear"
            run_overlapper
            run_trim_galore
            run_assembly_spades
            report_quast
            report_cegma
             busco_v2_dbs="eukaryota_odb9"
            report_busco_v2
            report_blobtools
            report_multiqc
            ;;
        h)
            help_message
            ;;
        \?)
            echo -e "Unknown Option."
            help_message
            ;;
    esac
done

if [ ! "$output_dir" ]
then
    echo "No output directory set. Please indicate -o" >&2
    exit 1
fi

exit 0
