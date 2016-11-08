#!/bin/bash

THREADS=18
SEQ_DIR="/storage/single_cells/completed"
DIR=$1

# match mito names to sequencing runs
# in an associative array
declare -A samples
samples=(["Katablepharid1_mt"]="katablepharids/Sample_11B_35C" ["Katablepharid3_mt"]="katablepharids/Sample_5F_35A" ["Katablepharid4_mt_contigs"]="katablepharids/Sample_6E_35B" ["Rhizaria14_mt.fa"]="rhizaria_97/Sample_8H_35C" ["Rhizaria15_mt"]="rhizaria_picozoa/Sample_10H_34B" ["Rhizaria6_mt"]="rhizaria_97/Sample_10C_34A" ["Rhizaria17_mt"]="rhizaria_picozoa/Sample_9G_34B" ["Rhizaria19_mt"]="cercozoa_90s/Sample_10G_35B" ["Rhizaria23_mt"]="cercozoa_90s/Sample_8G_35B" ["Rhizaria24_mt"]="cercozoa_other/Sample_2E_35B" ["Rhizaria27_mt"]="cryothecomonas/Sample_10D_34B" ["Rhizaria28_mt"]="cryothecomonas/Sample_4G_35C" ["Rhizaria29_mt"]="cryothecomonas/Sample_4H_35A" ["Rhizaria30_mt"]="cryothecomonas/Sample_5H_34A" ["Rhizaria31_mt"]="cryothecomonas/Sample_6D_35A" ["Rhizaria6_mt"]="rhizaria_97/Sample_10C_34A" ["Rhizaria8_mt.fa"]="rhizaria_97/Sample_3E_35B" ["Stramenopile10_FUNGAL_mt"]="mast_3/Sample_11F_35B" ["Stramenopile10_mt"]="mast_3/Sample_11F_35B" ["Stramenopile11_mt"]="mast_3/Sample_4D_34A" ["Stramenopile12_mt"]="mast_3/Sample_7A_34A" ["Stramenopile14_mt"]="mast_other/Sample_10G_35C" ["Stramenopile15_mt"]="mast_other/Sample_10H_35B" ["Stramenopile16_mt"]="mast_other/Sample_2E_34B" ["Stramenopile1_mt"]="stramenopiles/Sample_10B_35B" ["Stramenopile2_mt"]="stramenopiles/Sample_4H_35C" ["Stramenopile3_mt"]="stramenopiles/Sample_5B_35A" ["Stramenopile4_mt"]="stramenopiles_Sample_8H_35B" ["Stramenopile6_mt"]="mast_12B/Sample_9G_34A" ["Stramenopile8_mt"]="mast_1C/Sample_11F_35A" ["Telonemid10_mt"]="telonemida_99_100/Sample_7E_35A" ["Telonemid11_mt"]="telonemida_other/Sample_3E_35A" ["Telonemid12_mt"]="telonemida_other/Sample_3G_36B" ["Telonemid13_mt"]="telonemida_other/Sample_7C_34B" ["Telonemid1_mt"]="telonemida_94_95/Sample_10B_35C" ["Telonemid2_mt"]="telonemida_94_95/Sample_5D_35A" ["Telonemid5_mt"]="telonemida_94_95/Sample_8F_35C" ["Telonemid6_mt"]="telonemida_99_100/Sample_11H_34A")

MT_GENOMES=(*.fa)

for file in $(find ${DIR} -name '*.fa')
do
  current_sample=$(basename $file .fa)
  file=$(basename $file)

  echo "Working on ${current_sample}"

  cd $DIR\/$current_sample

  echo "Indexing Assembly $file"
  bwa index -a bwtsw $file | tee $file\_bwa.log

  # map reads to assembly with BWA MEM
  # we have to do this for all 5 sets of reads and then merge
  echo "Mapping Assembled reads to $file"
  bwa mem -t $THREADS $file \
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
  rm *.sam

  cd ../
done
