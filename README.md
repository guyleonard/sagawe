# Single Amplified Genome Assembly Workflow Example (SAG-AWE)

A suggested workflow for the assembly of [MDA](https://en.wikipedia.org/wiki/Multiple_displacement_amplification) sequenced Single Amplified Genomes from Illumina Hi-/Mi-Seq paired-end libraries.

## SAG Assembly Workflow Diagram
![SAGA Workflow](https://cdn.rawgit.com/guyleonard/single_cell_workflow/master/images/single_cell_workflow.svg)

## Example Usage / Help
Program parameters are order based; e.g. '-n' must come before '-a' and '-p' must come before '-t' etc...

Dependecnies are checked before running, please make sure they are in your PATH.

NB - There are several paths, relating to database locations, that must be manually set in the top of the run_single_cell_assemblies.sh script - these vary depending on your system set up.

    Single Amplified Genome Assembly Workflow
    Required Option:
      -o <output_dir>	Output Directory
    File Options:
      -f <forward.fastq>	Forward Reads
      -r <reverse.fastq>	Reverse Reads
    Optional Parameters (ordered):
      -n 	Read Normalisation
      -S 	Use Scaffolds Instead of Contigs
      -a 	Run All Options Below (p{bbmerge}tsqcb{eukaryota_odb9}Bm)
    Workflow Parameters:
      -p <pear|bbmerge>	Overlap Reads
      -t 	Trim Overlapped Reads
      -s 	Assemble Trimmed Reads
    Reports:
      -q 	Run QUAST
      -c 	Run CEGMA
      -b <db1,db2,...>	Run BUSCO v2
      -B 	Run BlobTools
      -m 	Run MultiQC
    Legacy (soon to be deprecated):
      -l 	Run BUSCO v1 - legacy
    Example: run_single_cell_assemblies.sh -f r1.fastq -r r2.fastq -o output_dir -n -S -a

### Output
Standard output, some folders/files have been truncated. Folder structure is assumed...

    output/
    ├── normalised
    ├── overlapped
    │   ├── assembled.fastq.gz
    │   ├── discarded.fastq.gz
    │   ├── unassembled.forward.fastq.gz
    │   └── unassembled.reverse.fastq.gz
    ├── trimmed
    │   ├── assembled_trimmed.fq.gz
    │   ├── unassembled.forward_val_1.fq.gz
    │   ├── unassembled.reverse_val_2.fq.gz
    │   ├── unassembled.forward_unpaired_1.fq.gz
    │   ├── unassembled.reverse_unpaired_2.fq.gz
    │   └── fastqc
    ├── assembly
    │   ├── contigs.fasta
    │   └── scaffolds.fasta
    └── reports
        ├── blobtools
        │   ├── blast
        │   │   └── scaffolds_vs_nt_1e-10.megablast
        │   ├── images
        │   │   ├── png
        │   │   └── svg
        │   ├── mapping
        │   │   └── scaffolds_mapped_all_reads.bam
        │   └── table
        ├── busco
        │   ├── run_bacteria_odb9
        │   │   └── short_summary_bacteria_odb9.txt
        │   ├── run_eukaryota_odb9
        │   │   └── short_summary_eukaryota_odb9.txt
        │   ├── run_protists_ensembl
        │   │   └── short_summary_protists_ensembl.txt
        │   └── summaries
        │       └── busco_figure.png
        ├── cegma
        │   └── cegma.completeness_report
        ├── multiqc_report.html
        └── quast
            ├── report.html
            ├── report.pdf
            └── report.txt

## SAG Assembly Workflow
1. Read Normalisation - Optional!! [BBNORM](http://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbnorm-guide/)
2. Overlapping of Read Libraries [BBMERGE](http://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmerge-guide/) or [PEAR](http://sco.h-its.org/exelixis/web/software/pear/doc.html) but since it's not easily accesible anymore, I will phase it out.
3. Trimming and Adaptor Cleaning [Trim Galore!](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
  * Remove any Illumina Sequencing Adaptors, poly-A tails, sequence quality score <20, etc.
4. Assembly of Prepared Reads [SPAdes](http://bioinf.spbau.ru/en/spades)
  * We find that SPAdes gives good results in SC mode, you may also like to try [IDBA-UD](http://i.cs.hku.hk/~alse/hkubrg/projects/idba_ud/index.html) and [Velvet-SC](http://bix.ucsd.edu/projects/singlecell/) both of which can assemble SCs. Not implemented here.
5. Assembly Statistics [QUAST](http://bioinf.spbau.ru/quast)
  * e.g. N50, contig/scaffold lengths/quantity, etc.
6. Read Mapping [BWA](https://github.com/lh3/bwa)
  * This can help with coverage information etc, but we will be using it mostly for "blobology" (see step 7)
7. BLAST Report
  * Top hits to NCBIs nt database using megablast, beware false top-hits as with any BLAST search.
8. BLOBology - GC/Taxonomy Maps [blobtools](https://github.com/DRL/blobtools)
  * GC/Coverage plots with taxonomy information to look for contamination.
9. Genome 'Completeness' Tests [CEGMA](http://korflab.ucdavis.edu/datasets/cegma/) & [BUSCO v1](http://busco.ezlab.org/v1/) - V2.0 Coming Soon!
10. MultiQC - Aggregate results from bioinformatics analyses across many samples into a single report

## Other Thoughts on Assembly & Downstream/Other Analyses
1. ESOM?
2. Co-assembly?
  * [HyDA](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4876485/)
  * Concat seq libraries + use pipeline with normalisation?
2. Merging Scaffolds (Reconcillition instead of Co-assembly??)
  * [Metassembler](https://sourceforge.net/projects/metassembler/)
  * Comparative Analysis and Merging of Scaffold Assemblies - [CAMSA](https://cblab.org/camsa/)
  * Contig Integrator for Sequence Assembly - [CISA](http://sb.nhri.org.tw/CISA/en/CISA)
3. Protocol for fully automated Decontamination of Genomes - [ProDeGe](http://www.nature.com/ismej/journal/v10/n1/full/ismej2015100a.html)
4. [CheckM](https://ecogenomics.github.io/CheckM/) - CheckM provides a set of tools for assessing the quality of genomes recovered from isolates, single cells, or metagenomes.
  * Looks like a nice tool, huge selection of options though and bac/arch oriented but can do ~euk.
  
## Citation
This work was initially completed for [in prep - paper here] for which the original scripts are available as a pre-release with the below DOI.
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.192677.svg)](https://doi.org/10.5281/zenodo.192677)

The repository and scripts however, have changed quite significantly since the initial release, any further citations should be to the below DOI.
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.438690.svg)](https://doi.org/10.5281/zenodo.438690)
