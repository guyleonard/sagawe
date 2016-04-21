#!/bin/bash
# Guy Leonard MMXVI

##
# This is a suggested workflow, it works on our servers/data...YMMV
##

##
# User Variables
##

# Number of processor cores
THREADS=8

## NCBI Databases
# NCBI 'nt' Database Location
NCBI_NT=/storage/ncbi/nt/nt
# NCBI Taxonomy Location
NCBI_TAX=/storage/ncbi/taxdump

## CEGMA Environment Variables
# CEGMA DIR
export CEGMA=/home/cs02gl/programs/CEGMA_v2
export PERL5LIB=$PERL5LIB:/home/cs02gl/programs/CEGMA_v2/lib
export WISECONFIGDIR=/usr/share/wise/
# Ammend PATH for CEGMA bin
export PATH=$PATH:/home/cs02gl/programs/CEGMA_v2/bin

## BUSCO Environment Variables
# BUSCO Lineage Location
BUSCO_DB=/storage/databases/BUSCO/eukaryota
# Augustus Config Path
export AUGUSTUS_CONFIG_PATH=~/programs/augustus-3.0.2/config/

## Locations of binaries if not in path
#SPADES
SPADES=/home/cs02gl/programs/SPAdes-3.7.0-Linux/bin/
#QUAST
QUAST=/home/cs02gl/programs/quast-4.0
#CEGMA
CEGMA_DIR=/home/cs02gl/programs/CEGMA_v2/bin
#BUSCO
BUSCO=/home/cs02gl/programs/BUSCO_v1.1b1
#BLOBTOOLS
BLOBTOOLS=/home/cs02gl/programs/blobtools

## Check programs are executable
command -v pigz >/dev/null 2>&1 || { echo "I require pigz but it's not installed.  Aborting." >&2; exit 1;}
command -v trim_galore >/dev/null 2>&1 || { echo "I require Trim Galore! but it's not installed.  Aborting." >&2; exit 1;}
command -v pear >/dev/null 2>&1 || { echo "I require PEAR but it's not installed.  Aborting." >&2; exit 1;}
command -v $SPADES/spades.py >/dev/null 2>&1 || { echo "I require SPAdes but it's not installed.  Aborting." >&2; exit 1;}
command -v $QUAST/quast.py >/dev/null 2>&1 || { echo "I require QUAST but it's not installed.  Aborting." >&2; exit 1;}
command -v $CEGMA_DIR/cegma >/dev/null 2>&1 || { echo "I require CEGMA but it's not installed.  Aborting." >&2; exit 1;}
command -v $BUSCO/BUSCO_v1.1b1.py >/dev/null 2>&1 || { echo "I require BUSCO but it's not installed.  Aborting." >&2; exit 1;}
command -v bwa >/dev/null 2>&1 || { echo "I require bwa but it's not installed.  Aborting." >&2; exit 1;}
command -v samtools1.3 >/dev/null 2>&1 || { echo "I require Samtools 1.3 but it's not installed.  Aborting." >&2; exit 1;}
command -v blastn >/dev/null 2>&1 || { echo "I require BLASTn but it's not installed.  Aborting." >&2; exit 1;}
command -v $BLOBTOOLS/blobtools >/dev/null 2>&1 || { echo "I require BLOBTOOLS but it's not installed.  Aborting." >&2; exit 1;}
command -v multiqc >/dev/null 2>&1 || { echo "I require MultiQC but it's not installed.  Aborting." >&2; exit 1;}

## Try not to change code below here...
# Working Directory
WD=`pwd`
echo "$WD"

# Get filenames for current Single Cell Library
# Locations of FASTQs = Sample_**_***/raw_illumina_reads/
for DIRS in */ ; do
	echo "Working in $DIRS"
	cd $DIRS/raw_illumina_reads

	# GZIP FASTQs
	# saving space down the line, all other files will be gzipped
	echo "gzipping *.fastq files"
	time pigz -9 -R *.fastq

	# Get all fastq.gz files
	FASTQ=(*.fastq.gz)

        # Run PEAR
        # default settings
        # output: pear_overlap
        mkdir -p PEAR
        cd PEAR
        echo "Running PEAR"
        time pear -f $WD/$DIRS/raw_illumina_reads/${FASTQ[0]} \
        -r $WD/$DIRS/raw_illumina_reads/${FASTQ[1]} \
        -o pear_overlap -j $THREADS | tee pear.log

        # Lets GZIP these too!
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
	time trim_galore -q 20 --fastqc --gzip --length 150 \
	$WD/$DIRS/raw_illumina_reads/PEAR/pear_overlap.assembled.fastq.gz

        echo "Running Trimming on Untrimmed Un-assembled Reads"
        time trim_galore -q 20 --fastqc --gzip --length 150 --paired --retain_unpaired \
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
	time $SPADES/spades.py --sc --careful -t $THREADS \
	--s1 $WD/$DIRS/raw_illumina_reads/PEAR/pear_overlap.assembled_trimmed.fq.gz \
        --s2 $WD/$DIRS/raw_illumina_reads/PEAR/pear_overlap.unassembled.forward_unpaired_1.fq.gz \
        --s3 $WD/$DIRS/raw_illumina_reads/PEAR/pear_overlap.unassembled.reverse_unpaired_2.fq.gz \
	--pe1-1 $WD/$DIRS/raw_illumina_reads/PEAR/pear_overlap.unassembled.forward_val_1.fq.gz \
	--pe1-2 $WD/$DIRS/raw_illumina_reads//PEAR/pear_overlap.unassembled.reverse_val_2.fq.gz \
	-o overlapped_and_paired | tee spades.log
	cd ../

	# Run QUAST
	# eukaryote mode
	# glimmer protein predictions
	mkdir -p QUAST
	cd QUAST
	echo "Running QUAST"
	time python $QUAST/quast.py -o quast_reports -t $THREADS \
	--min-contig 100 -f --eukaryote --scaffolds \
	--glimmer $WD/$DIRS/raw_illumina_reads/SPADES/overlapped_and_paired/scaffolds.fasta | tee quast.log
	cd ../

	# Run CEGMA
	# Genome mode
	echo "Running CEGMA"
	mkdir -p CEGMA
	cd CEGMA
	time $CEGMA_DIR/cegma -T $THREADS -g $WD/$DIRS/raw_illumina_reads/SPADES/overlapped_and_paired/scaffolds.fasta -o cegma
	cd ../

	# Run BUSCO
	echo "Running BUSCO"
	mkdir -p BUSCO
	cd BUSCO
	python3 $BUSCO/BUSCO_v1.1b1.py \
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
	time bwa index -a bwtsw $WD/$DIRS/raw_illumina_reads/SPADES/overlapped_and_paired/scaffolds.fasta | tee bwa.log

	# map reads to assembly with BWA MEM
	# we will have to do this for all 5 sets of reads and then merge
	echo "Mapping Assembled reads to Assembly"
	time bwa mem -t $THREADS $WD/$DIRS/raw_illumina_reads/SPADES/overlapped_and_paired/scaffolds.fasta \
	$WD/$DIRS/raw_illumina_reads/PEAR/pear_overlap.assembled_trimmed.fq.gz \
	> $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_assembled_reads.sam | tee -a bwa.log

        echo "Mapping Un-assembled reads to Assembly"
        time bwa mem -t $THREADS $WD/$DIRS/raw_illumina_reads/SPADES/overlapped_and_paired/scaffolds.fasta \
	$WD/$DIRS/raw_illumina_reads/PEAR/pear_overlap.unassembled.forward_unpaired_1.fq.gz \
	$WD/$DIRS/raw_illumina_reads/PEAR/pear_overlap.unassembled.reverse_unpaired_2.fq.gz \
        > $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_unassembled_reads.sam | tee -a bwa.log

        echo "Mapping Un-paired reads to Assembly"
        time bwa mem -t $THREADS $WD/$DIRS/raw_illumina_reads/SPADES/overlapped_and_paired/scaffolds.fasta \
	$WD/$DIRS/raw_illumina_reads/PEAR/pear_overlap.unassembled.forward_val_1.fq.gz \
	$WD/$DIRS/raw_illumina_reads/PEAR/pear_overlap.unassembled.reverse_val_2.fq.gz \
        > $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_unpaired_reads.sam | tee -a bwa.log

        # sort and convert sam to bam with SAMTOOLS
        #echo "Sorting Assembled SAM File and Converting to BAM"
        time samtools1.3 sort -@ $THREADS -o $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_assembled_reads.bam \
        $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_assembled_reads.sam | tee -a samtools.log

        # sort and convert sam to bam with SAMTOOLS
        #echo "Sorting Un-assembled SAM File and Converting to BAM"
        time samtools1.3 sort -@ $THREADS -o $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_unassembled_reads.bam \
        $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_unassembled_reads.sam | tee -a samtools.log

        # sort and convert sam to bam with SAMTOOLS
        #echo "Sorting Un-paired SAM File and Converting to BAM"
        time samtools1.3 sort -@ $THREADS -o $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_unpaired_reads.bam \
        $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_unpaired_reads.sam | tee -a samtools.log

	# Merge SAM files
	echo "Merging 3 BAM files"
	time samtools1.3 merge -n -@ $THREADS $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_all_reads.bam *.bam

	echo "Indexing Bam"
	time samtools1.3 index $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_all_reads.bam | tee -a samtools.log

	if [ ! -f $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_all_reads.bam.bai]
	then
		echo -e "[ERROR]\t[$DIRS]: No index file was created for your BAM file. !?" >> $WD/$DIRS/raw_illumina_reads/errors.txt
		# blobtools create will crash without this file, so we might as well move on to the next library...
		break
	fi

	# delete sam file - save some disk space, we have the bam now
	rm *.sam
	cd ../

	# run blast against NCBI 'nt'
	mkdir -p BLAST
	cd BLAST
	echo "Running BLAST"
	time blastn -task megablast \
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
	time $BLOBTOOLS/blobtools create -i $WD/$DIRS/raw_illumina_reads/SPADES/overlapped_and_paired/scaffolds.fasta \
	--nodes $NCBI_TAX/nodes.dmp --names $NCBI_TAX/names.dmp \
	-t $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/BLAST/scaffolds_vs_nt_1e-10.megablast \
	-b $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_all_reads.bam \
	-o scaffolds_mapped_reads_nt_1e-10_megablast_blobtools | tee -a $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/blobtools.log

	if  [ ! -f $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/scaffolds_mapped_reads_nt_1e-10_megablast_blobtools.BlobDB.json]
	then
		echo -e "[ERROR]\t[$DIRS]: Missing blobtools JSON, no tables or figures produced." >> $WD/$DIRS/raw_illumina_reads/errors.txt
	else
		# run blobtools view - table output
		# Standard Output - Phylum
		echo "Running BlobTools View"
		time $BLOBTOOLS/blobtools view -i $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/scaffolds_mapped_reads_nt_1e-10_megablast_blobtools.BlobDB.json \
		--out $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/scaffolds_mapped_reads_nt_1e-10_megablast_blobtools_phylum_table.csv | tee -a blobtools.log
		# Other Output - Species
		time $BLOBTOOLS/blobtools view -i $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/scaffolds_mapped_reads_nt_1e-10_megablast_blobtools.BlobDB.json \
		--out $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/scaffolds_mapped_reads_nt_1e-10_megablast_blobtools_superkingdom_table.csv \
		--rank superkingdom | tee -a blobtools.log

		# run blobtools plot - image output
		# Standard Output - Phylum, 7 Taxa
		echo "Running BlobTools Plots - Standard + SVG"
		time $BLOBTOOLS/blobtools plot -i $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/scaffolds_mapped_reads_nt_1e-10_megablast_blobtools.BlobDB.json
		time $BLOBTOOLS/blobtools plot -i $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/scaffolds_mapped_reads_nt_1e-10_megablast_blobtools.BlobDB.json \
		--format svg | tee -a blobtools.log

		# Other Output - Species, 15 Taxa
		echo "Running BlobTools Plots - SuperKingdom + SVG"
		time $BLOBTOOLS/blobtools plot -i $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/scaffolds_mapped_reads_nt_1e-10_megablast_blobtools.BlobDB.json \
		-r superkingdom
		time $BLOBTOOLS/blobtools plot -i $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/scaffolds_mapped_reads_nt_1e-10_megablast_blobtools.BlobDB.json \
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
done
