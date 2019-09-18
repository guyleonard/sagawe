# Single Amplified Genome - Assembly Workflow Example (SAG-AWE)

This script is intended as a suggested workflow for the assembly and analysis of [MDA](https://en.wikipedia.org/wiki/Multiple_displacement_amplification) sequenced 'Single Amplified Genomes' derived from Illumina Hi-/Mi-Seq paired-end libraries.

There are very few options available to the user, and most underlying programs are set with their defaults. This is intended to produce quick rough and reproducible assemblies of 100s of SAGs at a time.

## SAG Assembly Workflow Example Diagram
Depending on the options you select, there are four paths to assembling your data.

* T (green) - Trimmed
* TN (blue) - Trimmed and Normalised
* TM (yellow) - Trimmed and Merged
* TNM (purple) - Trimmed, Normalised and Merged

We suggest using all three options -t, -n and -m to produce the 'best' assembly, however your library prep/design and sequencing results may work better with different options. You can run the workflow with different options in the same output directory if you wish to make comparisons. 

<p align="center">
<img src="https://github.com/guyleonard/sagawe/blob/master/images/SAGAWE.svg">
</p>

* Red arrows indicate data requirements from previous steps
* Dotted arrows indicate reports read by MultiQC
* Grey arrows indicate future additions.

## Example Usage / Help
```
Single Amplified Genome Assembly Workflow Example (SAG-AWE)
  Options are positional, i.e. they are run sequentially, e.g. -S must come before -q.
Input Options (required):
  -f <r1.fq|r1.fq.gz> Read Library Pair 1
  -r <r2.fq|r2.fq.gz> Read Library Pair 2
Output Options (required):
  -o <output_dir> Output Directory
Program Parameters:
  -t  Run Trim Galore!
  -n  Run Normalisation
  -m  Run Read Merging
  -s  Run Assembly
Optional Parameters:
  -C  Use to turn off Single-Cell mode in SPAdes 
  -S  Use scaffolds.fasta instead of contigs.fasta
Reports/Stats:
  -k  Run KAT Analysis (run to inform QUAST)
  -q  Run QUAST Analysis
  -c  Run CEGMA Analysis
  -b </path/to/db1,/path/to/db2,...>  Run BUSCO with Multiple Lineages
  -B </path/to/blast/db,/path/to/taxdump> Run Blobtools
  -g  Run GenomeScope (requires -k)
  -p  Run Smudge Plots (requires -k and Jellyfish)
  -Q  Run Qualimap Analysis
  -M  Run MultiQC Analysis

Example: bin/sagawe -f read1.fq.gz -r read2.fq.gz -o results -t -n -m -s
```
The *-C* mode allows you to turn off "single-cell" mode in SPAdes, this may be useful for some analyses.

## Install Dependencies
The following programs will need to be installed, and be accessible from your PATH. The script has been tested on Ubuntu Linux Xenial, however it should work in other -nix environments. You will also need to make sure that the relevant environment variables for BLAST, CEGMA & AUGUSTUS are set correctly. This script does not check for dependencies and may fail without warning. Trim Galore! and SPAdes are the minimal required toolset you will need to start.

You may like to ask your local friendly bioinformatician / sys-admin to install the following programs for you in exchange for some beer or chocolate! :p Or you may try to install them yourself with conda.

### Assembly
* Trim_Galore! (required)
  * cutadapt
  * FastQC
  * pigz
* bbmap
  * bbmerge
  * bbnorm
* SPAdes (required)

### Reporting
* QUAST v4 or v5
  * glimmer
  * NCBI BLAST+
* Blobtools v1
  * BWA
  * samtools
  * NCBI BLAST+
* BUSCO v3
  * AUGUSTUS
  * bamtools
  * HMMER
* KAT
* Genomescope
* Smudgeplots
  * Jellyfish
* MultiQC

### Legacy
* CEGMA

You may like to try and install many of the dependencies via 'conda'. YMMV, the order below seemed to work for me, however it will give you warnings whilst solving the environment but in the end it should figure it out...

    conda create --name sags python=2.7
    conda activate sags
    conda install -c bioconda multiqc jellyfish kat blast busco blobtools bwa quast spades trim-galore bbmap

You will also need:
* [BUSCO Lineage Datasets](https://busco.ezlab.org)
* NCBI 'nt' Database
  * wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.??.tar.gz
* NCBI 'taxdump' Database
  * wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
* NCBI 'taxdb' Database
  * wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
  * Environment variable BLASTDB should be exported to the correct path
    * export BLASTDB=/path/to/taxdb
* [GenomeScope](https://github.com/tbenavi1/genomescope2.0)
* [SmudgePlot](https://github.com/KamilSJaron/smudgeplot)

NB - CEGMA is no longer supported and does not have a conda install, it is included for legacy purposes only, please use BUSCO.

## Citation
This work was initially designed during the making of [in prep]. The original scripts are available as a pre-release with the below DOI [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.192677.svg)](https://doi.org/10.5281/zenodo.192677).

The script has been rewritten quite significantly since the initial release and so any further citations should also include DOI.
Coming Soon