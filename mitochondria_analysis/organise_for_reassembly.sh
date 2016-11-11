#!/bin/bash

for dir in *; do
  if [[ -d $dir ]]; then
    echo "Working in ${dir}"
    cd ${dir}

    mkdir -p raw_illumina_reads/PEAR

    cp $PWD/${dir}_unpaired_reads.fastq $PWD/raw_illumina_reads/PEAR/pear_overlap.unpaired.fq
    cp $PWD/${dir}_all_reads.fastq.1 $PWD/raw_illumina_reads/PEAR/pear_overlap.paired.forward.fq
    cp $PWD/${dir}_all_reads.fastq.2 $PWD/raw_illumina_reads/PEAR/pear_overlap.paired.reverse.fq

    for file in $PWD/raw_illumina_reads/PEAR/*.fq; do
      pigz -9 -R "$file"
    done

    cd ../
  fi
done
