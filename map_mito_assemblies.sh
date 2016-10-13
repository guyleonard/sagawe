#!/bin/bash

THREADS=18
SEQ_DIR="/storage/single_cells/completed"

declare -A samples
# telonemids and katablepharids
samples=(["Telonemid2_mt"]="telonemida_94_95/Sample_10B_35C" ["Telonemid5_mt"]="telonemida_94_95/Sample_8F_35C" ["Telonemid6_mt"]="telonemida_99_100/Sample_11H_34A" ["Telonemid10_mt"]="telonemida_99_100/Sample_7E_35A" ["Telonemid11_mt"]="telonemida_other/Sample_3E_35A" ["Telonemid12_mt"]="telonemida_other/Sample_3G_36B" ["Telonemid13_mt"]="telonemida_other/Sample_7C_34B" ["Katablepharid1_mt_contigs"]="katablepharids/Sample_11B_35C" ["Katablepharid2_mt_contigs"]="katablepharids/Sample_11H_35C" ["Katablepharid4_mt_contigs"]="katablepharids/Sample_6E_35B")

# rhizaria
samples=(["Rhizaria14_mt.fa"]="rhizaria_97/Sample_8H_35C" ["Rhizaria16_mt.fa"]="rhizaria_picozoa/Sample_7E_34B" ["Rhizaria19_mt.fa"]="rhizaria_cercozoa_90s/Sample_10G_35B" ["Rhizaria24_mt.fa"]="rhizaria_cercozoa_other/Sample_2E_35B" ["Rhizaria28_mt.fa"]="rhizaria_cryothecomonas/Sample_4G_35C" ["Rhizaria30_mt.fa"]="rhizaria_cryothecomonas/Sample_5H_34A" ["Rhizaria6_mt.fa"]="rhizaria_97/Sample_10C_34A" ["Rhizaria15_mt.fa"]="rhizaria_picozoa/Sample_10H_34B" ["Rhizaria17_mt.fa"]="rhizaria_picozoa/Sample_9G_34B" ["Rhizaria23_mt.fa"]="rhizaria_cercozoa_90s/Sample_8G_35B" ["Rhizaria27_mt.fa"]="rhizaria_cryothecomonas/Sample_10D_34B" ["Rhizaria29_mt.fa"]="rhizaria/cryothecomonas/Sample_4H_35A" ["Rhizaria31_mt.fa"]="rhizaria_cryothecomonas/Sample_6D_35A" ["Rhizaria8_mt.fa"]="rhizaria_97/Sample_3E_35B")

MT_GENOMES=(*.fa)

for file in *.fa
do
  current_sample=$(basename $file .fa)

  echo "Indexing Assembly $file"
  bwa index -a bwtsw $file | tee $file\_bwa.log
  
  # map reads to assembly with BWA MEM
  # we have to do this for all 5 sets of reads and then merge
  echo "Mapping Assembled reads to $file"
  bwa mem -t $THREADS  $file \
  $SEQ_DIR/${samples[$current_sample]}/raw_illumina_reads/PEAR/pear_overlap.assembled_trimmed.fq.gz \
  > $current_sample\_mapped_assembled_reads.sam | tee -a bwa.log
  
  echo "Mapping Un-assembled & Un-Paired reads to Assembly - Forward"
  bwa mem -t $THREADS $file \
  $SEQ_DIR/${samples[$current_sample]}/raw_illumina_reads/PEAR/pear_overlap.unassembled.forward_unpaired_1.fq.gz \
  > $current_sample\_mapped_unassembled_unpaired_forward_reads.sam | tee -a bwa.log
  
  echo "Mapping Un-assembled & Un-Paired reads to Assembly - Reverse"
  bwa mem -t $THREADS $file \
  $SEQ_DIR/${samples[$current_sample]}/raw_illumina_reads/PEAR/pear_overlap.unassembled.reverse_unpaired_2.fq.gz \
  > $current_sample\_mapped_unassembled_unpaired_reverse_reads.sam | tee -a bwa.log
  
  echo "Mapping Un-assembled but still Paired reads to Assembly"
  bwa mem -t $THREADS $file \
  $SEQ_DIR/${samples[$current_sample]}/raw_illumina_reads/PEAR/pear_overlap.unassembled.forward_val_1.fq.gz \
  $SEQ_DIR/${samples[$current_sample]}/raw_illumina_reads/PEAR/pear_overlap.unassembled.reverse_val_2.fq.gz \
  > $current_sample\_mapped_unassembled_paired_reads.sam | tee -a bwa.log
  
  # sort and convert sam to bam with SAMTOOLS
  echo "Sorting Assembled SAM File and Converting to BAM"
  samtools sort -@ $THREADS -o $current_sample\_mapped_assembled_reads.bam \
  $current_sample\_mapped_assembled_reads.sam | tee -a samtools.log
  
  # sort and convert sam to bam with SAMTOOLS
  echo "Sorting Un-assembled & Un-Paired SAM Files and Converting to BAM - Forward"
  samtools sort -@ $THREADS -o $current_sample\_mapped_unassembled_unpaired_forward_reads.bam \
  $current_sample\_mapped_unassembled_unpaired_forward_reads.sam | tee -a samtools.log
  
  echo "Sorting Un-assembled & Un-Paired SAM Files and Converting to BAM - Reverse"
  samtools sort -@ $THREADS -o $current_sample\_mapped_unassembled_unpaired_reverse_reads.bam \
  $current_sample\_mapped_unassembled_unpaired_reverse_reads.sam | tee -a samtools.log
  
  # sort and convert sam to bam with SAMTOOLS
  echo "Sorting Un-assembled but still Paired SAM File and Converting to BAM"
  samtools sort -@ $THREADS -o $current_sample\_mapped_unassembled_paired_reads.bam \
  $current_sample\_mapped_unassembled_paired_reads.sam | tee -a samtools.log

  # Merge SAM files
  echo "Merging 4 BAM files"
  samtools merge -@ $THREADS -f $current_sample\_mapped_all_reads.bam \
  $current_sample\_mapped_assembled_reads.bam \
  $current_sample\_mapped_unassembled_paired_reads.bam \
  $current_sample\_mapped_unassembled_unpaired_forward_reads.bam \
  $current_sample\_mapped_unassembled_unpaired_reverse_reads.bam | tee -a samtools.log
  
  echo "Indexing Bam"
  samtools index $current_sample\_mapped_all_reads.bam | tee -a samtools.log
  
  # delete sam file - save some disk space, we have the bam now
  #rm *.sam

done