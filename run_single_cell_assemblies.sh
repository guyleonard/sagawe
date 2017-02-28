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
BUSCO_DB="/home/ubuntu/busco/eukaryota"
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
    mkdir -p "$normalised_dir"

    echo "Running BBNorm"
    bbnorm.sh in1="$READ1" in2="$READ2" \
    out1="$normalised_dir/$READ1{}" \
    out2="$normalised_dir/$READ2{}" \
    outt="$normalised_dir/excluded_reads.fastq.gz" \
    hist="$normalised_dir/input_kmer_depth.hist" \
    histout="$normalised_dir/output_kmer_depth.hist" \
    threads="$THREADS"
}

## Function to control which overlapper is used
function run_overlapper () {
    if $overlap_option eq "bbmerge" ; then
        run_bbmerge
    elif $overlap_option eq "pear" ; then
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

    echo "Running PEAR Assembler"
    pear -f "$READ1" -r "$READ2" \
    -o "$overlapped_dir/pear_overlap" -j "$THREADS" | tee "$overlapped_dir/pear.log"

    absolute_path="$( cd "$overlapped_dir" && pwd )"

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

    echo "Running BBMerge"
    bbmerge.sh in1="$READ1" in2="$READ2" \
    out="$absolute_path/assembled.fastq.gz" \
    outu1="$absolute_path/unassembled.forward.fastq.gz" \
    outu2="$absolute_path/unassembled.reverse.fastq.gz" \
    minoverlap=10 ziplevel=9
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
    fi
}

## Run SPAdes
# single cell mode - default kmers 21,33,55
# careful mode - runs mismatch corrector
# use all 5 sets of output reads
# 1 x assembled
# 2 x unassembled
# 2 x unpaired
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
        --s2 "$trimmed_dir/unassembled.forward_unpaired_1.fq.gz" \
        --s3 "$trimmed_dir/unassembled.reverse_unpaired_2.fq.gz" \
        --pe1-1 "$trimmed_dir/unassembled.forward_val_1.fq.gz" \
        --pe1-2 "$trimmed_dir/unassembled.reverse_val_2.fq.gz" \
        -o "$assembly_dir"
    fi
}

# Run QUAST
# eukaryote mode
# glimmer protein predictions
function report_quast () {
    assembly_dir="$output_dir/assembly"

    if [ ! -f "$assembly_dir/scaffolds.fasta" ] ; then
        echo -e "[ERROR]: SPAdes scaffolds cannot be found. Aborting."
        exit 1
    else
        quast_dir="$output_dir/reports/quast"
        mkdir -p "$quast_dir"
        echo "Running QUAST"
        quast.py -o "$quast_dir" -t "$THREADS" \
        --min-contig 100 -f --eukaryote --scaffolds \
        --glimmer "$assembly_dir/scaffolds.fasta" | tee "$quast_dir/quast.log"
    fi
}

# Run CEGMA
# Genome mode
function report_cegma () {
    assembly_dir="$output_dir/assembly"

    if [ ! -f "$assembly_dir/scaffolds.fasta" ] ; then
        echo -e "[ERROR]: SPAdes scaffolds cannot be found. Aborting."
        exit 1
    else
        cegma_dir="$output_dir/reports/cegma"
        mkdir -p "$cegma_dir"

        echo "Running CEGMA"
        cegma -T "$THREADS" -g "$assembly_dir/scaffolds.fasta" -o "$cegma_dir/cegma" | tee "$cegma_dir/cegma.log"

        # Tidy up, CEGMA's -o option doesn't seem to output to a dir!?
        absolute_path="$( cd "$output_dir" && pwd )"
        mv "$absolute_path/cegma."* "$absolute_path/reports/cegma"
    fi
}

# Run BUSCO
function report_busco () {
    assembly_dir="$output_dir/assembly"

    if [ ! -f "$assembly_dir/scaffolds.fasta" ] ; then
        echo -e "[ERROR]: SPAdes scaffolds cannot be found. Aborting."
        exit 1
    else
        busco_dir="$output_dir/reports/busco"
        mkdir -p "$busco_dir"

        absolute_path="$( cd "$output_dir" && pwd )"

        BUSCO_v1.22.py \
        -g "$assembly_dir/scaffolds.fasta" \
        -c "$THREADS" -l "$BUSCO_DB" -t "$absolute_path/reports/busco" -o "$absolute_path/reports/busco" -f | tee "$busco_dir/busco.log"

        # Tidy up busco
        mv "$absolute_path/reports/busco."* "$absolute_path/reports/busco"
    fi
}

function report_multiqc () {
    multiqc -o "$output_dir/reports" -z "$output_dir"
}

function blobtools_bwa () {
    assembly_dir="$output_dir/assembly"

    if [ ! -f "$assembly_dir/scaffolds.fasta" ] ; then
        echo -e "[ERROR]: SPAdes scaffolds cannot be found. Aborting."
        exit 1
    else
        blobtools_dir="$output_dir/reports/blobtools"
        mkdir -p "$blobtools_dir"

        blobtools_map="$blobtools_dir/mapping"
        mkdir -p "$blobtools_map"

        absolute_path="$( cd "$output_dir" && pwd )"

        trimmed_dir="$absolute_dir/trimmed"

        ln -s "$absolute_path/assembly/scaffolds.fasta" "$absolute_path/reports/blobtools/mapping/scaffolds.fasta"

        # index assembly (scaffolds.fa) with BWA
        echo "Indexing Assembly"
        bwa index -a bwtsw "$blobtools_map/scaffolds.fasta" | tee "$blobtools_map/bwa.log"

        # map reads to assembly with BWA MEM
        # we will have to do this for all 5 sets of reads and then merge
        echo "Mapping Assembled reads to Assembly"
        bwa mem -t "$THREADS" "$blobtools_map/scaffolds.fasta" \
        "$trimmed_dir/assembled_trimmed.fq.gz" \
        > "$blobtools_map/scaffolds_mapped_assembled_reads.sam" | tee -a "$blobtools_map/bwa.log"

        echo "Mapping Un-assembled & Un-Paired reads to Assembly - Forward"
        bwa mem -t "$THREADS" "$blobtools_map/scaffolds.fasta" \
        "$trimmed_dir/unassembled.forward_unpaired_1.fq.gz" \
        > "$blobtools_map/scaffolds_mapped_unassembled_unpaired_forward_reads.sam" | tee -a "$blobtools_map/bwa.log"

        echo "Mapping Un-assembled & Un-Paired reads to Assembly - Reverse"
        bwa mem -t "$THREADS" "$blobtools_map/scaffolds.fasta" \
        "$trimmed_dir/unassembled.reverse_unpaired_2.fq.gz" \
        > "$blobtools_map/scaffolds_mapped_unassembled_unpaired_reverse_reads.sam" | tee -a "$blobtools_map/bwa.log"

        echo "Mapping Un-assembled but still Paired reads to Assembly"
        bwa mem -t "$THREADS" "$blobtools_map/scaffolds.fasta" \
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
    -query "$blobtools_map/scaffolds.fasta" \
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
    blobtools create -i "$blobtools_map/scaffolds.fasta" \
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
    if [ ! -d "$BUSCO_DB" ] ; then
        echo "[ERROR]: Missing BUSCO Lineage Directory. Is your path correct?"
        echo "$BUSCO_DB"
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
    echo -e "Single Amplified Genome Assembly Pipeline"
    echo -e "Basic Usage:"
    echo -e "Required Parameters:"
    echo -e "  -f <forward.fastq>"
    echo -e "  -r <reverse.fastq>"
    echo -e "  -o <./output_dir>"
    echo -e "Pipeline Parameters:"
    echo -e "  -a 	Run All Options Below (ptsqcbBm)"
    echo -e "  -p <pear|bbmerge>	Overlap Reads"
    echo -e "  -t 	Trim Overlapped Reads"
    echo -e "  -s   Assemble Trimmed Reads"
    echo -e "  -n   (use|perform) Read Normalisation"
    echo -e "Reports:"
    echo -e "  -q 	Run QUAST"
    echo -e "  -c 	Run CEGMA"
    echo -e "  -b 	Run BUSCO"
    echo -e "  -B 	Run BlobTools"
    echo -e "  -m 	Run MultiQC"
    echo -e "Example: run_single_cell_assemblies.sh -f r1.fastq -r r2.fastq -o output_dir -a -n"
    exit 1
}

##############################
## Initial Pipeline Actions ##
##############################

# Set Cores
THREADS=$(cores)

# Check that we have the required programs
exes=(pigz clumpify.sh trim_galore pear spades.py quast.py cegma BUSCO_v1.22.py bwa samtools blastn blobtools multiqc)
for program in ${exes[@]} ; do
    check_exe "$program"
done

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

######################
## Pipeline Options ##
######################
while getopts f:r:o:np:tsqcbBmah FLAG; do
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
            ;;
        n)
            run_normalisation
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
        b)
            report_busco
            ;;
        B)
            report_blobtools
            ;;
        m)
            report_multiqc
            ;;
        a)
            overlap_option=$OPTARG
            run_overlapper
            run_trim_galore
            run_assembly_spades
            report_quast
            report_cegma
            report_busco
            report_blobtools
            report_multiqc
            ;;
        h)
            help_message
            ;;
        \?)
            echo -e "Option not allowed."
            help_message
            ;;
    esac
done

exit 0
