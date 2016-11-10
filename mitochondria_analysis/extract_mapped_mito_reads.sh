#!/bin/bash

for dir in *; do
    if [[ -d $dir ]]; then
       echo "Working in ${dir}"
       cd ${dir}

       # remove old/previous files
       echo "Removing old .fastq"
       #rm *.fastq
       #rm *.fq

       # remove other old junk files
       #echo "Removing junk files"
       #rm "${dir}_mapped_all_reads.bam"
       #rm "${dir}_mapped_all_reads.bam.bai"

	# samtools index the bam files
	echo "Indexing *.bam files"
	samtools index "${dir}_mapped_assembled_reads.bam"
	samtools index "${dir}_mapped_unassembled_paired_reads.bam"
	samtools index "${dir}_mapped_unassembled_unpaired_forward_reads.bam"
	samtools index "${dir}_mapped_unassembled_unpaired_reverse_reads.bam"

	# get the scaffolds from the fasta file
        # and extract the reads from the bame
	echo "Extract Reads for each Scaffold"
	grep ">" ${dir}.fa | while read -r line ; do

		# by the power of magic pattern substitution!
		scaffold=${line/>/}
		node_array=(${line//_/ })
		node="${node_array[1]}_${node_array[2]}"

		echo "Extracting reads that match to ${node}"
		samtools view -bh ${dir}_mapped_assembled_reads.bam ${scaffold} > ${node}_mapped_assembled_reads.bam
		samtools view -bh ${dir}_mapped_unassembled_paired_reads.bam ${scaffold} > ${node}_mapped_unassembled_paired_reads.bam
		samtools view -bh ${dir}_mapped_unassembled_unpaired_forward_reads.bam ${scaffold} > ${node}_mapped_unassembled_unpaired_forward_reads.bam
		samtools view -bh ${dir}_mapped_unassembled_unpaired_reverse_reads.bam ${scaffold} > ${node}_mapped_unassembled_unpaired_reverse_reads.bam

		echo "Converting extracted bams to fastq"
		bamtools convert -format fastq -in ${node}_mapped_assembled_reads.bam > ${node}_mapped_assembled_reads.fq
		bamtools convert -format fastq -in ${node}_mapped_unassembled_unpaired_forward_reads.bam > ${node}_mapped_unassembled_unpaired_forward_reads.fq
		bamtools convert -format fastq -in ${node}_mapped_unassembled_unpaired_reverse_reads.bam > ${node}_mapped_unassembled_unpaired_reverse_reads.fq
		bamtools convert -format fastq -in ${node}_mapped_unassembled_paired_reads.bam > ${node}_mapped_unassembled_paired_unordered_interleaved_reads.fq
	done

	# concatenate all .fq reads to respective libraries
        echo "concat .fq files to .fastq"
	cat *assembled_reads.fq > ${dir}_assembled_reads.fastq
	cat *unassembled_paired_unordered_interleaved_reads.fq > ${dir}_unassembled_paired_unordered_interleaved_reads.fastq
	cat *unassembled_unpaired_forward_reads.fq > ${dir}_unassembled_unpaired_forward_reads.fastq
	cat *unassembled_unpaired_reverse_reads.fq > ${dir}_unassembled_unpaired_reverse_reads.fastq

        echo "Sorting interleaved fastq"
	# use khmer script to split out orphaned pairs
	split-paired-reads.py -0 orphans ${dir}_unassembled_paired_unordered_interleaved_reads.fastq
	# append the orphans to one of the single reads files
	cat orphans >> ${dir}_unassembled_unpaired_reverse_reads.fastq
	#rm orphans
	# re-inter-leave/lace the properly paired reads
	interleave-reads.py ${dir}_unassembled_paired_unordered_interleaved_reads.fastq.1 ${dir}_unassembled_paired_unordered_interleaved_reads.fastq.2 -o ${dir}_unassembled_paired_ordered_interleaved_reads.fastq

	#rm ${dir}_unassembled_paired_unordered_interleaved_reads.fastq
	#rm ${dir}_unassembled_paired_unordered_interleaved_reads.fastq.1
	#rm ${dir}_unassembled_paired_unordered_interleaved_reads.fastq.2

	mkdir -p node_fq
	mv *.fq node_fq

	mkdir -p node_bam
	mv NODE*.bam* node_bam

	cd ../
    fi
done
