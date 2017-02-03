#!/usr/bin/env bash
# Guy Leonard MMXVI

##
# This is a suggested workflow, this will work on a new Amazon AMI of Ubuntu Xenial
# after using the bundles install_dependancies.sh script.
##

##
# User Defined Variables
##

# Number of Processor Cores / 2
THREADS=($(cores) / 2)

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

# Check that we have the required programs
exes=('pigz' 'clumpify' 'trim_galore' 'pear' 'spades.py' 'quast.py' 'cegma' 'BUSCO_v1.2.py' 'bwa' 'samtools' 'blastn' 'blobtools' 'multiqc')
for program in "${exes[@]}" ; do
  check_exe "$program"
done

while getopts f:r:h FLAG; do
    case $FLAG in
        f)
	    READ1=$OPTARG
	    ;;
	r)
	    READ2=$OPTARG
	    ;;
	h)
	    help
	    ;;
	\?) #unrecognized option - show help
         echo -e \\n"Option -$OPTARG not allowed."
         help
         ;;
    esac
done


## Try not to change code below here...
# Working Directory
WD=$(pwd)
echo "$WD"

# Get filenames for current Single Cell Library
# Locations of FASTQs = Sample_**_***/raw_illumina_reads/
for DIRS in */ ; do
	echo "Working in $DIRS"
	current_dir=$WD/$DIRS/raw_illumina_reads
	
	# GZIP FASTQs
	# saving space down the line, all other files will be gzipped
	echo "gzipping *.fastq files"
	time pigz -9 -R ./*.fastq

	# Get all fastq.gz files
	FASTQ=(*.fastq.gz)

        # Run PEAR
        mkdir -p "$current_dir/PEAR"
	run_pear

        # Lets GZIP these!
        echo "gzipping fastq files"
        pigz -9 -R $WD/$DIRS/raw_illumina_reads/PEAR/*.fastq

	# Run Trim Galore!
	# We need to do two sets of trimming:
	# one on the assembled reads
	# one on the unassembled paired reads.
	# minimum length of 150
	# minimum quality of Q20
	# run FASTQC on trimmed
	# GZIP output
	echo "Running Trimming on Untrimmed Assembled Reads"
	trim_galore -q 20 --fastqc --gzip --length 150 \
	$WD/$DIRS/raw_illumina_reads/PEAR/pear_overlap.assembled.fastq.gz

        echo "Running Trimming on Untrimmed Un-assembled Reads"
        trim_galore -q 20 --fastqc --gzip --length 150 --paired --retain_unpaired \
        $WD/$DIRS/raw_illumina_reads/PEAR/pear_overlap.unassembled.forward.fastq.gz pear_overlap.unassembled.reverse.fastq.gz

	cd ../

	# Run SPAdes
	# single cell mode - default kmers 21,33,55
	# careful mode - runs mismatch corrector
        # use all 5 sets of output reads
	# 1 x assembled
	# 2 x unassembled
	# 2 x unpaired
	mkdir -p SPADES
	cd SPADES
	echo "Running SPAdes"
	spades.py --sc --careful -t $THREADS \
	--s1 $WD/$DIRS/raw_illumina_reads/PEAR/pear_overlap.assembled_trimmed.fq.gz \
        --s2 $WD/$DIRS/raw_illumina_reads/PEAR/pear_overlap.unassembled.forward_unpaired_1.fq.gz \
        --s3 $WD/$DIRS/raw_illumina_reads/PEAR/pear_overlap.unassembled.reverse_unpaired_2.fq.gz \
	--pe1-1 $WD/$DIRS/raw_illumina_reads/PEAR/pear_overlap.unassembled.forward_val_1.fq.gz \
	--pe1-2 $WD/$DIRS/raw_illumina_reads//PEAR/pear_overlap.unassembled.reverse_val_2.fq.gz \
	-o overlapped_and_paired | tee spades.log
	cd ../

	# on occasion SPAdes, even though it is aware of the memory limits, will request more memory than is available
        # and then crash, we don't want the rest of the workflow to run through, and it would be nice to have an error message
	if [ ! -f $WD/$DIRS/raw_illumina_reads/SPADES/overlapped_and_paired/scaffolds.fasta ] ; then
		echo -e "[ERROR]\t[$DIRS]: SPAdes did not build scaffolds. This is possibly a memory error. This will need re-running" >> $WD/$DIRS/raw_illumina_reads/errors.txt
	else

	# Run QUAST
	# eukaryote mode
	# glimmer protein predictions
	mkdir -p QUAST
	cd QUAST
	echo "Running QUAST"
	quast.py -o quast_reports -t $THREADS \
	--min-contig 100 -f --eukaryote --scaffolds \
	--glimmer $WD/$DIRS/raw_illumina_reads/SPADES/overlapped_and_paired/scaffolds.fasta | tee quast.log
	cd ../

	# Run CEGMA
	# Genome mode
	echo "Running CEGMA"
	mkdir -p CEGMA
	cd CEGMA
	cegma -T $THREADS -g $WD/$DIRS/raw_illumina_reads/SPADES/overlapped_and_paired/scaffolds.fasta -o cegma
	cd ../

	# Run BUSCO
	echo "Running BUSCO"
	mkdir -p BUSCO
	cd BUSCO
	BUSCO_v1.2.py \
        -g $WD/$DIRS/raw_illumina_reads/SPADES/overlapped_and_paired/scaffolds.fasta \
	-c $THREADS -l $BUSCO_DB -o busco -f
	cd ../

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

## Program Functions ##

function run_pear () {
    echo "Running PEAR"
    pear -f $current_dir/PEAR/${FASTQ[0]} \
    -r $current_dir/PEAR/${FASTQ[1]} \
    -o pear_overlap -j $THREADS | tee pear.log
}

## Accessory Functions ##
function check_exe () {
    program=$1
    #https://techalicious.club/tutorials/validating-if-external-program-exists-linuxunix-based-bash-and-perl-scripts
    command -v $program >/dev/null 2>&1 || { echo "$program is required, but it is not in PATH.  Please install/ammend." >&2; exit 1;}
}

function cores () {
    cores=$(grep -c ^processor /proc/cpuinfo 2>/dev/null || sysctl -n hw.ncpu)
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
    echo -e "Example: run_single_cell_assemblies.sh -f r1.fastq -r r2.fastq"
    exit 1
}
