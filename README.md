# Single Amplified Genome - Assembly Workflow Example (SAG-AWE)

This script is intended as a suggested workflow for the assembly and analysis of [MDA](https://en.wikipedia.org/wiki/Multiple_displacement_amplification) sequenced 'Single Amplified Genomes' derived from Illumina Hi-/Mi-Seq paired-end libraries.

There are very few options available to the user, and most underlying programs are set with their defaults. This is intended to produce quick rough and reproducible assemblies of 100s of SAGs at a time.

## SAG Assembly Workflow Example Diagram
Depending on the options you select, there are four paths to assembling your data.

* T (green) - Trimmed
* TN (blue) - Trimmed and Normalised
* TM (yellow) - Trimmed and Merged
* TNM (blue+red=purple) - Trimmed, Normalised and Merged

We suggest using all three options -t, -n and -m to produce the 'best' assembly, however your library prep/design and sequencing results may work better with different options. You can run the workflow with different options in the same output directory if you wish to make comparisons. 

![SAGAWE](https://github.com/guyleonard/sagawe/blob/devel/images/SAGAWE.svg)

* Dotted lines indicate reports read by MultiQC
* Grey lines indicate future additions.

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
      -B </path/to/blast/db,/path/to/taxdump>  Run Blobtools with NCBI BLAST db
      -M  Run MultiQC
    Example: sag_awe -f read1.fq.gz -r read2.fq.gz -o results -t -n -m -s

## Install Dependencies
The following programs will need to be installed and be accessible from your PATH. Ihave only tested this on Ubuntu Linux, it should work in other -nix environments. You will also need to make sure that the environment variables for BLAST, CEGMA & AUGUSTUS are set correctly. This script does not check for them and may fail without warning.

Below is a brief list of the bioinformatic tools that you will need, some have their own dependencies. You may need to ask your friendly bioinformatician / sys admin to install them for you in exchange for some beer or chocolate! :p Or try to install them yourself with conda.

### Assembly
* Trim_Galore!
  * cutadapt
  * FastQC
  * pigz
* bbtools
  * bbmerge
  * bbnorm
* SPAdes

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
* Preseq
* MultiQC

### Legacy
* CEGMA

You may like to try and install many of the dependencies via 'conda'. YMMV.

    conda install -c bioconda trim-galore # installs Trim_Galore!, cutadapt, FastQC
    conda install -c agbiome bbtools # installs bbtools, samtools
    conda install -c bioconda spades # installs spades
    conda install -c bioconda quast # installs quast, blast, glimmer
    conda install -c bioconda bwa # installs bwa
    conda install -c bioconda blast # installs blast
    conda install -c bioconda blobtools # installs blobtools, samtools
    conda install -c bioconda busco # installs augustus, bamtools, blast, busco, hmmer
    conda install -c bioconda preseq # installs preseq
    conda install -c bioconda multiqc # installs multiqc

You will also need:
* [BUSCO Lineage Datasets](https://busco.ezlab.org)
* [NCBI 'nt' Database](ftp://ftp.ncbi.nlm.nih.gov/blast/db/)
* [NCBI 'taxdump' Database](ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/)
* [NCBI 'taxdb' Database](ftp://ftp.ncbi.nlm.nih.gov/blast/db/) - environment variable BLASTDB should be exported
  * export BLASTDB=/path/to/taxdb

NB - CEGMA is no longer supported and does not have a conda install and is included for legacy, please use BUSCO.

## Citation
This work was initially designed during the making of [in prep]. The original scripts are available as a pre-release with the below DOI [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.192677.svg)](https://doi.org/10.5281/zenodo.192677).

The script has been rewritten quite significantly since the initial release and so any further citations should also include DOI.
Coming Soon