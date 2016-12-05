# Single Cell Genome Assembly Workflow

This is a suggested workflow for the assembly of [MDA](https://en.wikipedia.org/wiki/Multiple_displacement_amplification) sequenced Single Cells from the Illumina MiSeq platform. These are sometimes also described as Single-cell Amplified Genomes (SAGs).

## Initial Paper Citation
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.192677.svg)](https://doi.org/10.5281/zenodo.192677)

## SAG Assembly Workflow Diagram
![SAGA Workflow](https://github.com/guyleonard/single_cell_workflow/blob/master/single_cell_workflow.png)

## SAG Assembly Workflow Explanation
This description is mostly program agnostic, however, we have placed the programs that we found to be most useful during our analyses and projects in brackets. You are welcome to adapt the workflow to inlude your favourite<sup>TM</sup> programs. 

1. Overlapping of Read Libraries<sup>[1](#footnote1)</sup> [PEAR](http://sco.h-its.org/exelixis/web/software/pear/doc.html)
  * We know our reads are 250bp PE with sequence overlaps, you might want to remove/adapt this step if yours are not.
2. Trimming and Adaptor Sequencing [Trim Galore!](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
  * Remove any Illumina Sequencing Adaptors, poly-A tails, sequence quality score <20, etc
3. Assembly of Reads [SPAdes](http://bioinf.spbau.ru/en/spades)
  * We find that SPAdes gives good results in SC mode, you may also like to try [IDBA-UD](http://i.cs.hku.hk/~alse/hkubrg/projects/idba_ud/index.html) and [Velvet-SC](http://bix.ucsd.edu/projects/singlecell/) both of which can assemble SCs.
4. Assembly Statistics [QUAST](http://bioinf.spbau.ru/quast)
  * e.g. N50, contig/scaffold lengths/quantity, etc.
5. Read Mapping [BWA](https://github.com/lh3/bwa)
  * This can help with coverage information etc, but we will be using it mostly for "blobology" (see step 7)
6. BLAST Report
  * Top hits to NCBI's 'nt' database using 'megablast', beware false hits as with any BLAST search.
7. 'BLOBology'- GC/Taxonomy Maps [blobtools](https://github.com/DRL/blobtools)
  * GC/Coverage plots with taxonomy information to look for contamination.
8. Genome 'Completeness' Tests CEGMA & BUSCO
9. MultiQC - Aggregate results from bioinformatics analyses across many samples into a single report

### Running
The directory structure is assumed (read: adapted to our local situation) to be:

* WorkDir/Sample_1/raw_illumina_reads/*.fastq
* ...
* WorkDir/Sample_n/raw_illumina_reads/*.fastq

The script will iterate over each folder one by one...

        run_single_cell_assemblies.sh | tee output_log.txt


## Gene Prediction Workflow

A suggested workflow for predicting genes from your assembly.

![SAGA Workflow](https://github.com/guyleonard/single_cell_workflow/blob/master/gene_prediction.png)

## Other Thoughts on Assembly & Downstream/Other Analyses

1. ESOM?
2. Contig Integrator for Sequence Assembly - [CISA](http://sb.nhri.org.tw/CISA/en/CISA)
3. Protocol for fully automated Decontamination of Genomes - [ProDeGe](http://www.nature.com/ismej/journal/v10/n1/full/ismej2015100a.html)
4. [CheckM](https://ecogenomics.github.io/CheckM/) - CheckM provides a set of tools for assessing the quality of genomes recovered from isolates, single cells, or metagenomes.
  * Looks like a nice tool, huge selection of options though and bac/arch oriented but can do ~euk.

## Footnotes
<a name="footnote1">1</a>: Originally I had steps 1 and 2 in the reverse order (trimming and then overlapping), on some read libraries this caused issues (I think where paired reads would become unordered and so overlapping would not run), however I don't think this is the problem with the order. Reads should be overlapped first, prior to trimming, as we should end up with a set of longer reads - due to the better quality scores over-riding the lower qualities within the overlapped areas, where this would not have happened if the lower quality reads were already trimmed.
 
