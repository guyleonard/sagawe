#!/bin/bash
## Guy Leonard MMXVI
## All programs are assumed to be in PATH, please make sure this is the case ;)

# Number of Cores to Use
THREADS=8

## Change these to your severs locations
# NCBI 'nt' Database Location
NCBI_NT=/storage/ncbi/nt/nt
# NCBI Taxonomy
NCBI_TAX=/storage/ncbi/taxdump

WD=`pwd`
echo "$WD"

# Get filenames for current Single Cell Library
# Locations of FASTQs = Sample_**_***/raw_illumina_reads/
for DIRS in */ ; do
	echo "Working in $DIRS"

	cd $DIRS/raw_illumina_reads

	# GZIP FASTQs
	# saving space down the line, all other files will be gzipped
	#echo "gzipping *.fastq files"
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
	cd ../

	# Run SPAdes
	# single cell mode - default kmers 21,33,55
	# careful - runs mismatch corrector
	mkdir -p SPADES
	cd SPADES
	echo "Running SPAdes"
	time spades.py --sc --careful -t $THREADS \
	--s1 $WD/$DIRS/raw_illumina_reads/PEAR/pear_overlap.assembled.fastq \
	--pe1-1 $WD/$DIRS/raw_illumina_reads/PEAR/pear_overlap.unassembled.forward.fastq \
	--pe1-2 $WD/$DIRS/raw_illumina_reads//PEAR/pear_overlap.unassembled.reverse.fastq \
	--pe2-1 ../${FILENAME[0]} \
	--pe2-2 ../${FILENAME[1]} \
	-o overlapped_and_paired | tee spades.log
	cd ../

	# Run QUAST
	# eukaryote mode
	# glimmer protein predictions
	mkdir -p QUAST
	cd QUAST
	echo "Running QUAST"
	time python quast.py -o quast_reports -t $THREADS \
	--min-contig 100 -f --eukaryote \
	--glimmer ../SPADES/overlapped_and_paired/contigs.fasta | tee quast.log
	cd ../

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
	time bwa mem -t $THREADS $WD/$DIRS/raw_illumina_reads/SPADES/overlapped_and_paired/scaffolds.fasta $WD/$DIRS/raw_illumina_reads/${FILENAME[0]} $WD/$DIRS/raw_illumina_reads/${FILENAME[1]} > $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_reads.sam | tee -a bwa.log

	# sort and convert sam to bam with SAMTOOLS
	echo "Sorting Sam File"
	time samtools1.3 sort -@ $THREADS -o $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_reads.bam $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_reads.sam | tee -a samtools.log

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

	cd ../../../
	echo "`pwd`"
	echo "Complete Run, Next or Finish."
done
