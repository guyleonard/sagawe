#!/bin/bash
# Guy Leonard MMXVI

##
# This is a suggested workflow, this will work on a new Amazon AMI of Ubuntu Xenial
# after using the bundles install_dependancies.sh script.
##

##
# User Variables
##

# Number of processor cores
THREADS=8

## NCBI Databases
# NCBI 'nt' Database Location and name (no extension)
NCBI_NT=/storage/ncbi/nt/nt
# NCBI Taxonomy Location
NCBI_TAX=/storage/ncbi/taxdump

## CEGMA Environment Variables
# CEGMA DIR
export CEGMA=/home/cs02gl/build/CEGMA_v2.5
export PERL5LIB=$PERL5LIB:/home/cs02gl/build/CEGMA_v2.5/lib
export WISECONFIGDIR=/usr/share/wise/
# Ammend PATH for CEGMA bin
export PATH=$PATH:/home/cs02gl/build/CEGMA_v2.5/bin

## BUSCO Environment Variables
# BUSCO Lineage Location
BUSCO_DB=/storage/databases/BUSCO/eukaryota
# Augustus Config Path
export AUGUSTUS_CONFIG_PATH=/home/cs02gl/build/augustus-3.2.2/config

## Check programs are executable
command -v pigz >/dev/null 2>&1 || { echo "I require pigz but it's not installed.  Aborting." >&2; exit 1;}
command -v trim_galore >/dev/null 2>&1 || { echo "I require Trim Galore! but it's not installed.  Aborting." >&2; exit 1;}
command -v pear >/dev/null 2>&1 || { echo "I require PEAR but it's not installed.  Aborting." >&2; exit 1;}
command -v spades.py >/dev/null 2>&1 || { echo "I require SPAdes but it's not installed.  Aborting." >&2; exit 1;}
command -v quast.py >/dev/null 2>&1 || { echo "I require QUAST but it's not installed.  Aborting." >&2; exit 1;}
command -v cegma >/dev/null 2>&1 || { echo "I require CEGMA but it's not installed.  Aborting." >&2; exit 1;}
command -v BUSCO_v1.2.py >/dev/null 2>&1 || { echo "I require BUSCO but it's not installed.  Aborting." >&2; exit 1;}
command -v bwa >/dev/null 2>&1 || { echo "I require bwa but it's not installed.  Aborting." >&2; exit 1;}
command -v samtools >/dev/null 2>&1 || { echo "I require Samtools 1.3 but it's not installed.  Aborting." >&2; exit 1;}
command -v blastn >/dev/null 2>&1 || { echo "I require BLASTn but it's not installed.  Aborting." >&2; exit 1;}
command -v blobtools >/dev/null 2>&1 || { echo "I require BLOBTOOLS but it's not installed.  Aborting." >&2; exit 1;}
command -v multiqc >/dev/null 2>&1 || { echo "I require MultiQC but it's not installed.  Aborting." >&2; exit 1;}

## Try not to change code below here...
# Working Directory
WD=`pwd`
echo "$WD"

# Get filenames for current Single Cell Library
# Locations of FASTQs = Sample_**_***/raw_illumina_reads/
for DIRS in */ ; do
	echo "Working in $DIRS"
	cd $WD/$DIRS/raw_illumina_reads

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
	spades.py --only-assembler -t $THREADS \
	--s1 $WD/$DIRS/raw_illumina_reads/PEAR/pear_overlap.unpaired.fq.gz \
        --pe1-1 $WD/$DIRS/raw_illumina_reads/PEAR/pear_overlap.paired.forward.fq.gz \
        --pe1-2 $WD/$DIRS/raw_illumina_reads/PEAR/pear_overlap.paired.reverse.fq.gz \
	-o $WD/$DIRS/raw_illumina_reads/SPADES/overlapped_and_paired | tee $WD/$DIRS/raw_illumina_reads/SPADES/spades.log
        #--pe1-1 $WD/$DIRS/raw_illumina_reads/PEAR/pear_overlap.unassembled.forward_val_1.fq.gz \
        #--pe1-2 $WD/$DIRS/raw_illumina_reads//PEAR/pear_overlap.unassembled.reverse_val_2.fq.gz \
	cd ../

	# on occasion SPAdes, even though it is aware of the memory limits, will request more memory than is available
        # and then crash, we don't want the rest of the workflow to run through, and it would be nice to have an error message
	if [ ! -f $WD/$DIRS/raw_illumina_reads/SPADES/overlapped_and_paired/scaffolds.fasta ]
	    then
		echo -e "[ERROR]\t[$DIRS]: SPAdes did not build scaffolds. This is possibly a memory error. This will need re-running" >> $WD/$DIRS/raw_illumina_reads/errors.txt
		cd ../
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
	echo "Mapping reads to Assembly"
	bwa mem -t $THREADS $WD/$DIRS/raw_illumina_reads/SPADES/overlapped_and_paired/scaffolds.fasta \
	$WD/$DIRS/raw_illumina_reads/PEAR/pear_overlap.unpaired.fq.gz \
	> $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_unpaired_reads.sam | tee -a bwa.log

        bwa mem -t $THREADS $WD/$DIRS/raw_illumina_reads/SPADES/overlapped_and_paired/scaffolds.fasta \
        $WD/$DIRS/raw_illumina_reads/PEAR/pear_overlap.paired.forward.fq.gz \
        > $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_paired_forward_reads.sam | tee -a bwa.log

        bwa mem -t $THREADS $WD/$DIRS/raw_illumina_reads/SPADES/overlapped_and_paired/scaffolds.fasta \
        $WD/$DIRS/raw_illumina_reads/PEAR/pear_overlap.paired.reverse.fq.gz \
        > $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_paired_reverse_reads.sam | tee -a bwa.log

        # sort and convert sam to bam with SAMTOOLS
        echo "Sorting SAM Files and Converting to BAM"
        samtools sort -@ $THREADS -o $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_unpaired_reads.bam \
        $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_unpaired_reads.sam | tee -a samtools.log

        samtools sort -@ $THREADS -o $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_paired_forward_reads.bam \
        $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_paired_forward_reads.sam | tee -a samtools.log

        samtools sort -@ $THREADS -o $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_paired_reverse_reads.bam \
        $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_paired_reverse_reads.sam | tee -a samtools.log

	# Merge SAM files
	echo "Merging 4 BAM files"
	samtools merge -@ $THREADS -f $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_all_reads.bam \
	$WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_unpaired_reads.bam \
	$WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_paired_forward_reads.bam \
	$WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_paired_reverse_reads.bam | tee -a samtools.log

	echo "Indexing Bam"
	samtools index $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_all_reads.bam | tee -a samtools.log

	if [ ! -f $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/MAPPING/scaffolds_mapped_all_reads.bam.bai ]
	then
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

	if  [ ! -f $WD/$DIRS/raw_illumina_reads/BLOBTOOLS/scaffolds_mapped_reads_nt_1e-10_megablast_blobtools.BlobDB.json ]
	  then
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
