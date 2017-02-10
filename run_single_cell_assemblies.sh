#!/usr/bin/env bash
# Guy Leonard MMXVI

##
# This is a suggested workflow, this will work on a new Amazon AMI of Ubuntu Xenial
# after using the bundles install_dependancies.sh script.
##

##
# User Defined Variables
##

# NCBI 'nt' Database Location and name (no extension)
NCBI_NT=/home/ubuntu/blast/nt/nt

# NCBI Taxonomy Location
NCBI_TAX=/home/ubuntu/blast/taxonomy

# CEGMA DIR
CEGMA_DIR=/home/ubuntu/single_cell_workflow/build/CEGMA_v2.5
export_cegma

# BUSCO Lineage Location
BUSCO_DB=/home/ubuntu/busco/eukaryota

# Augustus Config Path
export AUGUSTUS_CONFIG_PATH=/home/ubuntu/single_cell_workflow/build/augustus-3.2.2/config

while getopts f:r:o:pth FLAG; do
    case $FLAG in
        f)
	        READ1=$OPTARG
	        ;;
	    r)
	        READ2=$OPTARG
	        ;;
	    o)
            $current_dir=$OPTARG
            ;;
	    p)
			overlapped_dir=$current_dir/overlapped
            mkdir -p "$overlapped_dir"
            run_pear
            ;;
        t)
			trimmed_dir=$current_dir/trimmed
			mkdir -p "$trimmed_dir"
		    trim_galore
		    ;;
		s)
			assembly_dir=$current_dir/assembly
			mkdir -p "$assembly_dir"
			assembly_spades
	    h)
	        help
	        ;;
	    \?)
            echo -e \\n"Option -$OPTARG not allowed."
            help
            ;;
    esac
done


## Try not to change code below here...
# Working Directory
#WD=$(pwd)
#echo "$WD"

# Get filenames for current Single Cell Library
# Locations of FASTQs = Sample_**_***/raw_illumina_reads/
for DIRS in */ ; do
	#echo "Working in $DIRS"
	#current_dir=$WD/$DIRS/raw_illumina_reads
	
	# GZIP FASTQs
	# saving space down the line, all other files will be gzipped
	#echo "gzipping *.fastq files"
	#time pigz -9 -R ./*.fastq

	# Get all fastq.gz files
	#FASTQ=(*.fastq.gz)

    # Lets GZIP these!
    #echo "gzipping fastq files"
    #pigz -9 -R $WD/$DIRS/raw_illumina_reads/PEAR/*.fastq


	# Run BlobTools
	mkdir -p BLOBTOOLS
	cd BLOBTOOLS
	mkdir -p MAPPING
	cd MAPPING
	# index assembly (scaffolds.fa) with BWA
	echo "Indexing Assembly"
	bwa index -a bwtsw $WD/$DIRS/raw_illumina_reads/SPADES/overlapped_and_paired/scaffolds.fasta | tee bwa.log

	# map reads to assembly with BWA MEM
	# we will have to do this for all 5 sets of reads and then merge
	echo "Mapping Assembled reads to Assembly"
	bwa mem -t $THREADS $WD/$DIRS/raw_illumina_reads/SPADES/overlapped_and_paired/scaffolds.fasta \
	$WD/$DIRS/raw_illumina_reads/PEAR/pear_overlap.assembled_trimmed.fq.gz \
	> $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_assembled_reads.sam | tee -a bwa.log

        echo "Mapping Un-assembled & Un-Paired reads to Assembly - Forward"
        bwa mem -t $THREADS $WD/$DIRS/raw_illumina_reads/SPADES/overlapped_and_paired/scaffolds.fasta \
	$WD/$DIRS/raw_illumina_reads/PEAR/pear_overlap.unassembled.forward_unpaired_1.fq.gz \
        > $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_unassembled_unpaired_forward_reads.sam | tee -a bwa.log

	echo "Mapping Un-assembled & Un-Paired reads to Assembly - Reverse"
        bwa mem -t $THREADS $WD/$DIRS/raw_illumina_reads/SPADES/overlapped_and_paired/scaffolds.fasta \
	$WD/$DIRS/raw_illumina_reads/PEAR/pear_overlap.unassembled.reverse_unpaired_2.fq.gz \
        > $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_unassembled_unpaired_reverse_reads.sam | tee -a bwa.log

        echo "Mapping Un-assembled but still Paired reads to Assembly"
        bwa mem -t $THREADS $WD/$DIRS/raw_illumina_reads/SPADES/overlapped_and_paired/scaffolds.fasta \
	$WD/$DIRS/raw_illumina_reads/PEAR/pear_overlap.unassembled.forward_val_1.fq.gz \
	$WD/$DIRS/raw_illumina_reads/PEAR/pear_overlap.unassembled.reverse_val_2.fq.gz \
        > $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_unassembled_paired_reads.sam | tee -a bwa.log

        # sort and convert sam to bam with SAMTOOLS
        echo "Sorting Assembled SAM File and Converting to BAM"
        samtools sort -@ $THREADS -o $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_assembled_reads.bam \
        $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_assembled_reads.sam | tee -a samtools.log

        # sort and convert sam to bam with SAMTOOLS
        echo "Sorting Un-assembled & Un-Paired SAM Files and Converting to BAM - Forward"
        samtools sort -@ $THREADS -o $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_unassembled_unpaired_forward_reads.bam \
        $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_unassembled_unpaired_forward_reads.sam | tee -a samtools.log

        echo "Sorting Un-assembled & Un-Paired SAM Files and Converting to BAM - Reverse"
        samtools sort -@ $THREADS -o $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_unassembled_unpaired_reverse_reads.bam \
        $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_unassembled_unpaired_reverse_reads.sam | tee -a samtools.log

        # sort and convert sam to bam with SAMTOOLS
        echo "Sorting Un-assembled but still Paired SAM File and Converting to BAM"
        samtools sort -@ $THREADS -o $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_unassembled_paired_reads.bam \
        $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_unassembled_paired_reads.sam | tee -a samtools.log

	# Merge SAM files
	echo "Merging 4 BAM files"
	samtools merge -@ $THREADS -f $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_all_reads.bam \
	$WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_assembled_reads.bam \
	$WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_unassembled_paired_reads.bam \
	$WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_unassembled_unpaired_forward_reads.bam \
	$WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_unassembled_unpaired_reverse_reads.bam | tee -a samtools.log

	echo "Indexing Bam"
	samtools index $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_all_reads.bam | tee -a samtools.log

	if [ ! -f $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_all_reads.bam.bai ] ; then
		echo -e "[ERROR]\t[$DIRS]: No index file was created for your BAM file. !?" >> $WD/$DIRS/raw_illumina_reads/errors.txt
		# blobtools create will crash without this file, so we might as well move on to the next library...
		#break
		#continue
	fi

	# delete sam file - save some disk space, we have the bam now
	rm $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/*.sam
	cd ../

	# run blast against NCBI 'nt'
	mkdir -p BLAST
	cd BLAST
	echo "Running BLAST"
	blastn -task megablast \
	-query $WD/$DIRS/raw_illumina_reads/SPADES/overlapped_and_paired/scaffolds.fasta \
	-db $NCBI_NT \
	-evalue 1e-10 \
	-num_threads $THREADS \
	-outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
	-culling_limit 5 \
	-out $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/BLAST/scaffolds_vs_nt_1e-10.megablast | tee blast.log
	cd ../

	# run blobtools create
	echo "Running BlobTools CREATE - slow"
	cd $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/
	blobtools create -i $WD/$DIRS/raw_illumina_reads/SPADES/overlapped_and_paired/scaffolds.fasta \
	--nodes $NCBI_TAX/nodes.dmp --names $NCBI_TAX/names.dmp \
	-t $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/BLAST/scaffolds_vs_nt_1e-10.megablast \
	-b $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_all_reads.bam \
	-o scaffolds_mapped_reads_nt_1e-10_megablast_blobtools | tee -a $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/blobtools.log

	if  [ ! -f $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/scaffolds_mapped_reads_nt_1e-10_megablast_blobtools.BlobDB.json ] ; then
		echo -e "[ERROR]\t[$DIRS]: Missing blobtools JSON, no tables or figures produced." >> $WD/$DIRS/raw_illumina_reads/errors.txt
	else
		# run blobtools view - table output
		# Standard Output - Phylum
		echo "Running BlobTools View"
		blobtools view -i $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/scaffolds_mapped_reads_nt_1e-10_megablast_blobtools.BlobDB.json \
		--out $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/scaffolds_mapped_reads_nt_1e-10_megablast_blobtools_phylum_table.csv | tee -a blobtools.log
		# Other Output - Species
		blobtools view -i $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/scaffolds_mapped_reads_nt_1e-10_megablast_blobtools.BlobDB.json \
		--out $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/scaffolds_mapped_reads_nt_1e-10_megablast_blobtools_superkingdom_table.csv \
		--rank superkingdom | tee -a blobtools.log

		# run blobtools plot - image output
		# Standard Output - Phylum, 7 Taxa
		echo "Running BlobTools Plots - Standard + SVG"
		blobtools plot -i $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/scaffolds_mapped_reads_nt_1e-10_megablast_blobtools.BlobDB.json
		blobtools plot -i $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/scaffolds_mapped_reads_nt_1e-10_megablast_blobtools.BlobDB.json \
		--format svg | tee -a blobtools.log

		# Other Output - Species, 15 Taxa
		echo "Running BlobTools Plots - SuperKingdom + SVG"
		blobtools plot -i $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/scaffolds_mapped_reads_nt_1e-10_megablast_blobtools.BlobDB.json \
		-r superkingdom
		blobtools plot -i $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/scaffolds_mapped_reads_nt_1e-10_megablast_blobtools.BlobDB.json \
		-r superkingdom \
		--format svg | tee -a blobtools.log
	fi

	# Run MultiQC for some extra, nice stats reports on QC etc
	cd ../
        multiqc $WD/$DIRS/raw_illumina_reads/

	# Finish up.
	cd ../../
	echo "`pwd`"
	echo "Complete Run, Next or Finish."

	# end if from scaffolds.fasta check
	fi
done


#######################
## Program Functions ##
#######################

## Run Pear Overlapper
# default settings
## deprecate for bbmerge?
function run_pear () {
    echo "Running PEAR Assembler"
    pear -f $overlapped_dir/$READ1 \
    -r $overlapped_dir/$READ2 \
    -o pear_overlap -j $THREADS | tee $overlapped_dir/pear.log

    ln -s $trimmed_dir/pear_overlap.assembled.fastq.gz $trimmed_dir/assembled.fastq.gz
    ln -s $trimmed_dir/pear_overlap.unassembled.forward.fastq.gz $trimmed_dir/unassembled.forward.fastq.gz
    ln -s $trimmed_dir/pear_overlap.unassembled.reverse.fastq.gz $trimmed_dir/unassembled.reverse.fastq.gz
}

## Run Trim Galore!
# Two sets of trimming:
# Assembled/overlapped reads
# Unassembled paired reads.
# Minimum length of 150
# Minimum quality of Q20
# run FASTQC on trimmed
# GZIP output
function trim_galore () {
	echo "Running Trim Galore! on Untrimmed Assembled Reads"
	trim_galore -q 20 --fastqc --gzip --length 150 \
	$trimmed_dir/assembled.fastq.gz

	echo "Running Trim Galore! on Untrimmed Un-assembled Reads"
	trim_galore -q 20 --fastqc --gzip --length 150 --paired --retain_unpaired \
	$trimmed_dir/unassembled.forward.fastq.gz \
	$trimmed_dir/unassembled.reverse.fastq.gz
}

## Run SPAdes
# single cell mode - default kmers 21,33,55
# careful mode - runs mismatch corrector
# use all 5 sets of output reads
# 1 x assembled
# 2 x unassembled
# 2 x unpaired
function assembly_spades () {
	echo "Running SPAdes"
	spades.py --sc --careful -t $THREADS \
	--s1 $trimmed_dir/assembled_trimmed.fq.gz \
	--s2 $trimmed_dir/unassembled.forward_unpaired_1.fq.gz \
	--s3 $trimmed_dir/unassembled.reverse_unpaired_2.fq.gz \
	--pe1-1 $trimmed_dir/unassembled.forward_val_1.fq.gz \
	--pe1-2 $trimmed_dir/unassembled.reverse_val_2.fq.gz \
	-o overlapped_and_paired | tee $assembly_dir/spades_overlapped_and_paired.log

	# Sometimes SPAdes, even though it is aware of the memory limits, will request more memory than is available
    # and then crash, we don't want the rest of the workflow to run through, and it would be nice to have an error message
	#if [ ! -f $assembly_dir/overlapped_and_paired/scaffolds.fasta ] ; then
	#	echo -e "[ERROR]: SPAdes did not build scaffolds. This is possibly a memory error."
	#	exit 1
	#fi
}

# Run QUAST
# eukaryote mode
# glimmer protein predictions
function report_quast () {

	if [ ! -f $assembly_dir/overlapped_and_paired/scaffolds.fasta ] ; then
		echo -e "[ERROR]: SPAdes scaffolds cannot be found."
		exit 1
	else
		echo "Running QUAST"
		quast.py -o $report_dir/quast -t $THREADS \
		--min-contig 100 -f --eukaryote --scaffolds \
		--glimmer $assembly_dir/overlapped_and_paired/scaffolds.fasta | tee quast.log
	fi
}

# Run CEGMA
# Genome mode
function report_cegma () {
	if [ ! -f $assembly_dir/overlapped_and_paired/scaffolds.fasta ] ; then
		echo -e "[ERROR]: SPAdes scaffolds cannot be found."
		exit 1
	else
		echo "Running CEGMA"
		cegma -T $THREADS -g $assembly_dir/overlapped_and_paired/scaffolds.fasta -o $report_dir/cegma
	fi
}


# Run BUSCO
function report_busco () {
	if [ ! -f $assembly_dir/overlapped_and_paired/scaffolds.fasta ] ; then
		echo -e "[ERROR]: SPAdes scaffolds cannot be found."
		exit 1
	else
		BUSCO_v1.22.py \
	    -g $assembly_dir/overlapped_and_paired/scaffolds.fasta \
		-c $THREADS -l $BUSCO_DB -o $report_dir/busco -f
	fi
}

function report_multiqc () {

}

function blobtools_bwa () {

}

function blobtools_blast () {

}

function blobtools_create () {

}

function blobtools_table () {

}

function blobtools_image () {

}


#########################
## Accessory Functions ##
#########################

function check_exe () {
    program=$1
    #https://techalicious.club/tutorials/validating-if-external-program-exists-linuxunix-based-bash-and-perl-scripts
    command -v $program >/dev/null 2>&1 || { echo "$program is required, but it is not in PATH.  Please install/ammend." >&2; exit 1;}
}

function cores () {
    cores=$(grep -c ^processor /proc/cpuinfo 2>/dev/null || sysctl -n hw.ncpu)
    echo $(($cores / 2))
}

function export_cegma () {
    export CEGMA=$CEGMA_DIR
    export PATH=$PATH:$CEGMA_DIR/bin
    export PERL5LIB=$PERL5LIB:$CEGMA_DIR/lib
    export WISECONFIGDIR=/usr/share/wise/
}

function HELP {
    echo -e "Basic Usage:"
    echo -e "-f Read 1 FASTQ"
    echo -e "-r Read 2 FASTQ"
    echo -e "-o Output Directory"
    echo -e "Example: run_single_cell_assemblies.sh -f r1.fastq -r r2.fastq -o output_dir"
    exit 1
}


## Main Pipeline Actions ##

THREADS=$(cores)
echo "num_threads:$THREADS"

# Check that we have the required programs
exes=('pigz' 'clumpify.sh' 'trim_galore' 'pear' 'spades.py' 'quast.py' 'cegma' 'BUSCO_v1.22.py' 'bwa' 'samtools' 'blastn' 'blobtools' 'multiqc')
for program in "${exes[@]}" ; do
    check_exe "$program"
done
