#!/bin/bash
# Guy Leonard MMXVI

## All programs are assumed to be in PATH, please make sure this is the case ;)
# This is just a suggested workflow, it works on our servers...you will have
# to adapt it to your location.

## This script uses CEGMA - notoriously difficult to install and now unsupported
# You might not want to use it, you are safe to comment it out

## It also uses BUSCO
# I still don't trust BUSCO hence also using CEGMA
# If you see something like this:
# Error: Sequence file ./run_testing//augustus_proteins/64334.fas is empty or misformatted
# just ignore it, BUSCO is 'working' it's an error from hmmer and there are no gene predictions
# for the SCOs

# Number of Cores to Use
THREADS=8

## Change these to your sever's directory locations

## NCBI
# NCBI 'nt' Database Location
NCBI_NT=/storage/ncbi/nt/nt
# NCBI Taxonomy
NCBI_TAX=/storage/ncbi/taxdump

## CEGMA
# CEGMA DIR
export CEGMA=/home/cs02gl/programs/CEGMA_v2
export PERL5LIB=$PERL5LIB:/home/cs02gl/programs/CEGMA_v2/lib

## BUSCO
# BUSCO Lineage Location
BUSCO=/storage/databases/BUSCO/eukaryota
export AUGUSTUS_CONFIG_PATH=~/programs/augustus-3.0.2/config/

## PATH
# This needs to have the CEGMA BIN directory added to it
export PATH=$PATH:/home/cs02gl/programs/CEGMA_v2

## Do not change below here...

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

	# Run Trim Galore!
	# minimum length of 150
	# minimum quality of Q20
	# run FASTQC on trimmed
	# GZIP output
	echo "Running Trimming"
	time trim_galore -q 20 --fastqc --gzip --length 150 \
	--paired $WD/$DIRS/raw_illumina_reads/${FASTQ[0]} $WD/$DIRS/raw_illumina_reads/${FASTQ[1]}

        # Get all fq.gz files - these are the default names from Trim Galore!
	# Making it nice and easy to distinguish from our original .fastq inputs
        FILENAME=(*.fq.gz)

	# Run PEAR
	# default settings
	# output: pear_overlap
	mkdir -p PEAR
	cd PEAR
	echo "Running PEAR"
	time pear -f $WD/$DIRS/raw_illumina_reads/${FILENAME[0]} \
        -r $WD/$DIRS/raw_illumina_reads/${FILENAME[1]} \
        -o pear_overlap -j $THREADS | tee pear.log

	# Lets GZIP these too!
	pigz -9 -R *.fastq
	cd ../

	# Run SPAdes
	# single cell mode - default kmers 21,33,55
	# careful mode - runs mismatch corrector
	mkdir -p SPADES
	cd SPADES
	echo "Running SPAdes"
	time spades.py --sc --careful -t $THREADS \
	--s1 $WD/$DIRS/raw_illumina_reads/PEAR/pear_overlap.assembled.fastq.gz \
	--pe1-1 $WD/$DIRS/raw_illumina_reads/PEAR/pear_overlap.unassembled.forward.fastq.gz \
	--pe1-2 $WD/$DIRS/raw_illumina_reads//PEAR/pear_overlap.unassembled.reverse.fastq.gz \
	#--pe2-1 ../${FILENAME[0]} \
	#--pe2-2 ../${FILENAME[1]} \
	-o overlapped_and_paired | tee spades.log
	cd ../

	# Run QUAST
	# eukaryote mode
	# glimmer protein predictions
	mkdir -p QUAST
	cd QUAST
	echo "Running QUAST"
	time python quast.py -o quast_reports -t $THREADS \
	--min-contig 100 -f --eukaryote --scaffolds \
	--glimmer $WD/$DIRS/raw_illumina_reads/SPADES/overlapped_and_paired/scaffolds.fasta | tee quast.log
	cd ../

	# Run CEGMA
	mkdir -p CEGMA
	#time cegma -T 8 -g $WD/$DIRS/raw_illumina_reads/SPADES/overlapped_and_paired/scaffolds.fasta -o testing

	# Run BUSCO
	mkdir -p BUSCO
	#python3 ~/programs/BUSCO_v1.1b1/BUSCO_v1.1b1.py \
        #-g $WD/$DIRS/raw_illumina_reads/SPADES/overlapped_and_paired/scaffolds.fasta \
	#-c 8 -l /storage/databases/BUSCO/eukaryota/ -o testing

	# Run BlobTools
	mkdir -p BLOBTOOLS
	cd BLOBTOOLS
	mkdir -p MAPPING
	cd MAPPING
	# index assembly (scaffolds.fa) with BWA
	echo "Indexing Assembly"
	time bwa index -a bwtsw $WD/$DIRS/raw_illumina_reads/SPADES/overlapped_and_paired/scaffolds.fasta | tee bwa.log

	# map original reads to assembly with BWA MEM
	echo "Mapping reads to Assembly"
	time bwa mem -t $THREADS $WD/$DIRS/raw_illumina_reads/SPADES/overlapped_and_paired/scaffolds.fasta $WD/$DIRS/raw_illumina_reads/${FILENAME[0]} \
	$WD/$DIRS/raw_illumina_reads/${FILENAME[1]} > $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_reads.sam | tee -a bwa.log

	# sort and convert sam to bam with SAMTOOLS
	echo "Sorting Sam File"
	time samtools1.3 sort -@ $THREADS -o $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_reads.bam \
	$WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_reads.sam | tee -a samtools.log

	echo "Converting Sam to Bam"
	time samtools1.3 index $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_reads.bam | tee -a samtools.log

	# delete sam file
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
	time blobtools create -i $WD/$DIRS/raw_illumina_reads/SPADES/overlapped_and_paired/scaffolds.fasta \
	--nodes $NCBI_TAX/nodes.dmp --names $NCBI_TAX/names.dmp \
	-t $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/BLAST/scaffolds_vs_nt_1e-10.megablast \
	-b $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_reads.bam \
	-o scaffolds_mapped_reads_nt_1e-10_megablast_blobtools | tee -a $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/blobtools.log

	# run blobtools view - table output
	# Standard Output - Phylum
	echo "Running BlobTools View"
	time blobtools/blobtools view -i $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/scaffolds_mapped_reads_nt_1e-10_megablast_blobtools.BlobDB.json \
	--out $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/scaffolds_mapped_reads_nt_1e-10_megablast_blobtools_phylum_table.csv | tee -a blobtools.log
	# Other Output - Species
	time blobtools/blobtools view -i $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/scaffolds_mapped_reads_nt_1e-10_megablast_blobtools.BlobDB.json \
	--out $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/scaffolds_mapped_reads_nt_1e-10_megablast_blobtools_superkingdom_table.csv \
	--rank superkingdom | tee -a blobtools.log

	# run blobtools plot - image output
	# Standard Output - Phylum, 7 Taxa
	echo "Running BlobTools Plots - Standard + SVG"
	time blobtools/blobtools plot -i $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/scaffolds_mapped_reads_nt_1e-10_megablast_blobtools.BlobDB.json
	time blobtools/blobtools plot -i $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/scaffolds_mapped_reads_nt_1e-10_megablast_blobtools.BlobDB.json \
	--format svg | tee -a blobtools.log

	# Other Output - Species, 15 Taxa
	echo "Running BlobTools Plots - SuperKingdom + SVG"
	time blobtools/blobtools plot -i $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/scaffolds_mapped_reads_nt_1e-10_megablast_blobtools.BlobDB.json \
	-r superkingdom
	time blobtools/blobtools plot -i $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/scaffolds_mapped_reads_nt_1e-10_megablast_blobtools.BlobDB.json \
	-r superkingdom \
	--format svg | tee -a blobtools.log

	# Finish up.
	cd ../../../
	echo "`pwd`"
	echo "Complete Run, Next or Finish."
done
