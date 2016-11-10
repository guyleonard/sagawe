#!/bin/bash
THREADS=18
SEQ_DIR="/storage/single_cells/completed"
MAP_DIR="/raw_illumina_reads/BLOBTOOLS/MAPPING"

# match mito names to sequencing runs
# in an associative array
declare -A samples
samples=(["Katablepharid1_mt"]="katablepharids/Sample_11B_35C" ["Katablepharid3_mt"]="katablepharids/Sample_5F_35A" ["Katablepharid4_mt"]="katablepharids/Sample_6E_35B" ["Rhizaria14_mt"]="rhizaria_97/Sample_8H_35C" ["Rhizaria15_mt"]="rhizaria_picozoa/Sample_10H_34B" ["Rhizaria16_mt"]="rhizaria_97/Sample_7E_34B" ["Rhizaria17_mt"]="rhizaria_picozoa/Sample_9G_34B" ["Rhizaria19_mt"]="cercozoa_90s/Sample_10G_35B" ["Rhizaria23_mt"]="cercozoa_90s/Sample_8G_35B" ["Rhizaria24_mt"]="cercozoa_other/Sample_2E_35B" ["Rhizaria27_mt"]="cryothecomonas/Sample_10D_34B" ["Rhizaria28_mt"]="cryothecomonas/Sample_4G_35C" ["Rhizaria29_mt"]="cryothecomonas/Sample_4H_35A" ["Rhizaria30_mt"]="cryothecomonas/Sample_5H_34A" ["Rhizaria31_mt"]="cryothecomonas/Sample_6D_35A" ["Rhizaria6_mt"]="rhizaria_97/Sample_10C_34A" ["Rhizaria8_mt"]="rhizaria_97/Sample_3E_35B" ["Stramenopile10_FUNGAL_mt"]="mast_3/Sample_11F_35B" ["Stramenopile10_mt"]="mast_3/Sample_11F_35B" ["Stramenopile11_mt"]="mast_3/Sample_4D_34A" ["Stramenopile12_mt"]="mast_3/Sample_7A_34A" ["Stramenopile14_mt"]="mast_other/Sample_10G_35C" ["Stramenopile15_mt"]="mast_other/Sample_10H_35B" ["Stramenopile16_mt"]="mast_other/Sample_2E_34B" ["Stramenopile1_mt"]="stramenopiles/Sample_10B_35B" ["Stramenopile2_mt"]="stramenopiles/Sample_4H_35C" ["Stramenopile3_mt"]="stramenopiles/Sample_5B_35A" ["Stramenopile4_mt"]="stramenopiles/Sample_8H_35B" ["Stramenopile6_mt"]="mast_12B/Sample_9G_34A" ["Stramenopile8_mt"]="mast_1C/Sample_11F_35A" ["Telonemid10_mt"]="telonemida_99_100/Sample_7E_35A" ["Telonemid11_mt"]="telonemida_other/Sample_3E_35A" ["Telonemid12_mt"]="telonemida_other/Sample_3G_36B" ["Telonemid13_mt"]="telonemida_other/Sample_7C_34B" ["Telonemid1_mt"]="telonemida_94_95/Sample_10B_35C" ["Telonemid2_mt"]="telonemida_94_95/Sample_5D_35A" ["Telonemid5_mt"]="telonemida_94_95/Sample_8F_35C" ["Telonemid6_mt"]="telonemida_99_100/Sample_11H_34A")

for dir in *; do
  if [[ -d $dir ]]; then
    echo "Working in ${dir}"
    cd $dir

    # remove old/previous files
    echo "Removing old .fastq files"
    rm *.fastq

    # samtools sort and index the bam files
    echo "Sorting *.bam file from ${SEQ_DIR}/${samples[$dir]}/${MAP_DIR}/scaffolds_mapped_all_reads_sorted.bam"
    samtools sort -@ $THREADS -o ${SEQ_DIR}/${samples[$dir]}/${MAP_DIR}/scaffolds_mapped_all_reads_sorted.bam ${SEQ_DIR}/${samples[$dir]}/${MAP_DIR}/scaffolds_mapped_all_reads.bam
    samtools index $SEQ_DIR/${samples[$dir]}/${MAP_DIR}/scaffolds_mapped_all_reads_sorted.bam

    # get the scaffolds from the fasta file
    # and extract the reads from the bame
    echo "Extract Reads for each Scaffold"
    grep ">" ${dir}.fa | while read -r line ; do

      # by the power of magic pattern substitution!
      node_array=(${line//_/ })
      node="${node_array[1]}_${node_array[2]}"
      scaffold="${node_array[1]}_${node_array[2]}_${node_array[3]}_${node_array[4]}_${node_array[5]}_${node_array[6]}_${node_array[7]}_${node_array[8]}"

      echo "Extracting reads that match to ${node}"
      samtools view -b $SEQ_DIR/${samples[$dir]}/${MAP_DIR}/scaffolds_mapped_all_reads_sorted.bam ${scaffold} > ${node}_mapped_reads.bam

      echo "Converting extracted bams to fastq"
      bamtools convert -format fastq -in ${node}_mapped_reads.bam > ${node}_mapped_reads.fq
    done

    # concatenate all .fq reads to respective libraries
    echo "concat .fq files to .fastq"
    cat *_reads.fq > ${dir}_all_reads.fastq

    echo "Sorting interleaved fastq"
    # use khmer script to split out orphaned pairs
    split-paired-reads.py -0 ${dir}_unpaired_reads.fastq ${dir}_all_reads.fastq

    mkdir -p nodes_fq
    mv *.fq nodes_fq

    mkdir -p nodes_bam
    mv NODE*.bam* nodes_bam

    cd ../
    fi
done
