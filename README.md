# Single Amplified Genome Assembly Workflow Example (SAG-AWE)

A suggested workflow for the assembly of [MDA](https://en.wikipedia.org/wiki/Multiple_displacement_amplification) sequenced Single Amplified Genomes from Illumina Hi-/Mi-Seq paired-end libraries.

A list of dependencies appears in the list describing the workflow. These programs will need to be installed correctly and accessible on the command line. You will also need to make sure that the environment variables for BLAST, CEGMA, AUGUSTUS etc are correctly set. This script does not check for them and may fail without.

## SAG Assembly Workflow Example Diagram
![SAGAWE](https://github.com/guyleonard/sagawe/blob/devel/images/SAGAWE.svg)
  
  * Depending on the options you selected, there are four paths to assembling your data.
    * T (green) - Trimmed - all reads need adaptor trimming
    * TN (blue) - Trimmed and Normalised
    * TM (yellow) - Trimmed and Merged
    * TNM (blue+red=purple) - Trimmed, Normalised and Merged (default/suggested/preffered)
  * Dotted lines indicate reports read by FastQC
  * Grey line indicates future addition

## Example Usage / Help
    Single Amplified Genome Assembly Workflow Example (SAG-AWE)
      Options are positional, i.e. they are run sequentially, e.g. -S must come before -q, -t before -n, etc.
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

## Install Dependencies
You may like to try and install many of the dependencies via 'conda'. YMMV.

    conda install -c bioconda trim-galore # installs Trim_Galore!, cutadapt, FastQC
    conda install -c agbiome bbtools # installs bbtools, samtools
    conda install -c bioconda spades # installs spades
    conda install -c bioconda quast # installs quast, blast, glimmer
    conda install -c bioconda bwa # installs bwa
    conda install -c bioconda blast # installs blast
    conda install -c bioconda blobtools # installs blobtools, samtools
    conda install -c bioconda busco # installs augustus, bamtools, blast, busco, hmmer
    conda install -c bioconda multiqc # installs multqiqc

* You will also need:
  * [BUSCO Lineage Datasets](https://busco.ezlab.org)
  * [NCBI 'nt' Database](ftp://ftp.ncbi.nlm.nih.gov/blast/db/)
  * [NCBI 'taxdump' Database](ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/)
  * [NCBI 'taxdb' Database](ftp://ftp.ncbi.nlm.nih.gov/blast/db/) - environment variable BLASTDB
    * export BLASTDB=/path/to/taxdb
* NB - CEGMA is no longer supported and does not have a conda install.

## Citation
This work was initially completed for [in prep - paper here] for which the original scripts are available as a pre-release with the below DOI.
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.192677.svg)](https://doi.org/10.5281/zenodo.192677)

The repository and scripts however, have changed quite significantly since the initial release, any further citations should be to the below DOI.
Coming Soon
