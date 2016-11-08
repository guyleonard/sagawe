#!/bin/bash

for dir in *; do
    if [[ -d $dir ]]; then
       echo "${dir}"
       cd $dir

       # remove old/previous files
       echo "Removing .fastq"
       rm *.fastq
       rm *.fq

       # remove other old junk files
       echo "Removing junk files"
       rm "${dir}_mt_mapped_unassembled_paired_reads.qsort.bam"
       rm "${dir}_mt_mapped_all_reads.bam"
       rm "${dir}_mt_mapped_all_reads.bam.bai"

	# samtools index the bam files
	echo "Indexing *.bam files"
	samtools index "${dir}_mt_mapped_assembled_reads.bam"
	samtools index "${dir}_mt_mapped_unassembled_paired_reads.bam"
	samtools index "${dir}_mt_mapped_unassembled_unpaired_forward_reads.bam"
	samtools index "${dir}_mt_mapped_unassembled_unpaired_reverse_reads.bam"

	# get the scaffolds from the fasta file
        # and extract the reads from the bame
	echo "Extract Reads for each Scaffold"
	grep ">" ${dir}_mt.fa | while read -r line ; do

		# by the power of magic pattern substitution!
		scaffold=${line/>/}
		node_array=(${line//_/ })
		node="${node_array[1]}_${node_array[2]}"

		echo "Extracting reads that match to ${node}"
		samtools view -bh ${dir}_mt_mapped_assembled_reads.bam ${scaffold} > ${node}_mapped_assembled_reads.bam
		samtools view -bh ${dir}_mt_mapped_unassembled_paired_reads.bam ${scaffold} > ${node}_mapped_unassembled_paired_reads.bam
		samtools view -bh ${dir}_mt_mapped_unassembled_unpaired_forward_reads.bam ${scaffold} > ${node}_mapped_unassembled_unpaired_forward_reads.bam
		samtools view -bh ${dir}_mt_mapped_unassembled_unpaired_reverse_reads.bam ${scaffold} > ${node}_mapped_unassembled_unpaired_reverse_reads.bam

		echo "Converting extracted bams to fastq"
		bamtools convert -format fastq -in ${node}_mapped_assembled_reads.bam > ${node}_mapped_assembled_reads.fq
		bamtools convert -format fastq -in ${node}_mapped_unassembled_unpaired_forward_reads.bam > ${node}_mapped_unassembled_unpaired_forward_reads.fq
		bamtools convert -format fastq -in ${node}_mapped_unassembled_unpaired_reverse_reads.bam > ${node}_mapped_unassembled_unpaired_reverse_reads.fq

		# bedtools can't handle pairs that are out of order, so we will have to have an interleaved fastq
		echo "Converting paired read bam to interleaved fastq"
		bamtools convert -format fastq -in ${node}_mapped_unassembled_paired_reads.bam > ${node}_mapped_unassembled_paired_unordered_interleaved_reads.fq
	done

	# concatenate all .fq reads to respective libraries
        echo "concat .fq files to .fastq"
	cat *assembled_reads.fq > ${dir}_assembled_reads.fastq
	cat *unassembled_paired_unordered_interleaved_reads.fq > ${dir}_unassembled_paired_unordered_interleaved_reads.fastq
	cat *unassembled_unpaired_forward_reads.fq > ${dir}_unassembled_unpaired_forward_reads.fastq
	cat *unassembled_unpaired_reverse_reads.fq > ${dir}_unassembled_unpaired_reverse_reads.fastq

        echo "Sorting interleaved fastq"
        #cat ${dir}_unassembled_paired_unordered_interleaved_reads.fastq | paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' > ${dir}_unassembled_paired_interleaved_reads.fastq
        #rm ${dir}_unassembled_paired_unordered_interleaved_reads.fastq

	# the sort still sometimes gets out of order
	# use khmer split-paired-reads.py -0 erroneous_unpaired and dump them in to one of the unpaired files
	# then use the khmer interleave-reads.py to put the corrected paired reads together

	mkdir node_fq
	mv *.fq node_fq

	mkdir node_bam
	mv NODE*.bam node_bam

	cd ../
    fi
done
