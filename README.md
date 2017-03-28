# Single Amplified Genome Assembly Workflow

A suggested workflow for the assembly of [MDA](https://en.wikipedia.org/wiki/Multiple_displacement_amplification) sequenced Single Amplified Genomes from Illumina Hi-/Mi-Seq paired-end libraries.

## SAG Assembly Workflow Diagram
![SAGA Workflow](https://cdn.rawgit.com/guyleonard/single_cell_workflow/master/images/single_cell_workflow.svg)

## Example Usage / Help
Program parameters are order based; e.g. -n must come before -a and -p before -t etc...
NB - There are seven paths, relating to database locations, that must be set in the top of the run_single_cell_assemblies.sh script - these vary depending on your system set up.

    Single Amplified Genome Assembly Pipeline
    Basic Usage:
    Required Parameters:
      -f <forward.fastq>
      -r <reverse.fastq>
      -o <./output_dir>
    Pipeline Parameters:
      -n  Read Normalisation (optional)
      -S  Use scaffolds instead of contigs
      -a  Run All Options Below (ptsqcbBm)
      -p  <pear|bbmerge>  Overlap Reads (pear is default)
      -t  Trim Overlapped Reads
      -s  Assemble Trimmed Reads
    Reports:
      -q 	Run QUAST
      -c 	Run CEGMA
      -b 	Run BUSCO
      -B 	Run BlobTools
      -m 	Run MultiQC
    
    Example: run_single_cell_assemblies.sh -f r1.fastq -r r2.fastq -o output_dir -n -a

## SAG Assembly Workflow
1. Read Normalisation - Optional!! [BBNORM](http://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbnorm-guide/)
2. Overlapping of Read Libraries [PEAR](http://sco.h-its.org/exelixis/web/software/pear/doc.html) or [BBMERGE](http://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmerge-guide/)
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

## Gene Prediction Workflow

A suggested workflow for predicting genes from your assembly.

![SAGA Workflow](https://github.com/guyleonard/single_cell_workflow/blob/master/images/gene_prediction.png)

## Other Thoughts on Assembly & Downstream/Other Analyses

1. ESOM?
2. Contig Integrator for Sequence Assembly - [CISA](http://sb.nhri.org.tw/CISA/en/CISA)
3. Protocol for fully automated Decontamination of Genomes - [ProDeGe](http://www.nature.com/ismej/journal/v10/n1/full/ismej2015100a.html)
4. [CheckM](https://ecogenomics.github.io/CheckM/) - CheckM provides a set of tools for assessing the quality of genomes recovered from isolates, single cells, or metagenomes.
  * Looks like a nice tool, huge selection of options though and bac/arch oriented but can do ~euk.
  
## Initial Paper Citation
This work was initially started from [insert paper here] for which the original scripts are available as a release with the below DOI. The repository and scripts have subsequently changed quite significantly, although the workflow remains much the same.
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.192677.svg)](https://doi.org/10.5281/zenodo.192677)
