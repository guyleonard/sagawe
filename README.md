# Single Amplified Genome Assembly Workflow Example (SAG-AWE)

A suggested workflow for the assembly of [MDA](https://en.wikipedia.org/wiki/Multiple_displacement_amplification) sequenced Single Amplified Genomes from Illumina Hi-/Mi-Seq paired-end libraries.

A list of dependencies appears in the list describing the workflow. These programs will need to be installed correctly and accessible on the command line. You will also need to make sure that the environment variables for BLAST, CEGMA, AUGUSTUS etc are correctly set. This script does not check for them and may fail without.

## Example Usage / Help
Program parameters are order based; e.g. '-n' must come before '-a' and '-p' must come before '-t' etc...

    Single Amplified Genome Assembly Workflow Example (SAG-AWE)
      Options are positional, i.e. they are run sequentially, e.g. -S must come before -q.
    Input Options (required):
      -f <r1.fq|r1.fq.gz>  Read Library Pair 1
      -r <r2.fq|r2.fq.gz>  Read Library Pair 2
    Output Options (required):
      -o <output_dir>  Output Directory
    Program Parameters:
      -t  Run Trim Galore!
      -n  Run Normalisation (bbnorm)
      -m  Run Read Merging (bbmerge)
      -s  Run Assembly
    Optional Parameters:
      -S  Use scaffolds.fasta instead of contigs.fasta
    Reports/Stats:
      -q  Run Quast
      -c  Run CEGMA
      -b </path/to/db1,/path/to/db2,...>  Run BUSCO with Multiple Lineages
      -B </path/to/blast/db>  Run Blobtools with NCBI BLAST db
    Example: sag_awe -f read1.fq.gz -r read2.fq.gz -o results -t -n -m -s

## SAG Assembly Workflow Example
1. Trimming and Adaptor Cleaning [Trim Galore!](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
  * Remove any Illumina Sequencing Adaptors, poly-A tails, sequence quality score <20, etc.
2. Read Normalisation - Optional!! [BBNORM](http://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbnorm-guide/)
  * Normalise your data - especially useful for SAGs.
3. Overlap Read Libraries [BBMERGE](http://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmerge-guide/)
  * Construct longer 'reads' and/or increases quality of reads.
4. Assembly of Prepared Reads [SPAdes](http://bioinf.spbau.ru/en/spades)
  * Depending on the options you selected in Steps 1-4
  ** T - Trimmed - all reads need adaptor trimming
  ** TN - Trimmed and Normalised
  ** TM - Trimmed and Merged
  ** TNM - Trimmed, Normalised and Merged (default/suggested/preffered)
5. Assembly Statistics [QUAST](http://bioinf.spbau.ru/quast)
  * N50/L50, contig/scaffold lengths/quantity, etc.
6. Read Mapping [BWA](https://github.com/lh3/bwa)
  * This can help with coverage information etc, but we will be using it mostly for "blobology" (see step 8)
7. BLAST Report
  * Top hits to NCBI style database using megablast.
8. BLOBology - GC/Taxonomy Maps [blobtools](https://github.com/DRL/blobtools)
  * GC/Coverage plots with taxonomy information to look for contamination.
9. Genome 'Completeness' Tests [CEGMA](http://korflab.ucdavis.edu/datasets/cegma/) & [BUSCO v1](http://busco.ezlab.org/v3/)
10. MultiQC - Aggregate results from bioinformatics analyses across many samples into a single report

## Install Dependencies

You may like to try and install many of the dependencies via 'conda'

    conda install -c bioconda trim-galore

## Other Information

## Citation
This work was initially completed for [in prep - paper here] for which the original scripts are available as a pre-release with the below DOI.
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.192677.svg)](https://doi.org/10.5281/zenodo.192677)

The repository and scripts however, have changed quite significantly since the initial release, any further citations should be to the below DOI.
Coming Soon
