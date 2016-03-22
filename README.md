# Single Cell Assembly Workflow

This is a suggested workflow for the assembly of MDA Single Cell MiSeq Illumina libraries. It is contained within a bash script that will execute all the programs sequentially.
The directory structure is assumed (read: adapted to our local situation) to be:

* WD/Sample_1/raw_illumina_reads/*.fastq
* WD/Sample_2/raw_illumina_reads/*.fastq

The script will iterate over each folder one by one...

## Workflow in Steps
 * PIGZ - Parallel GZIP *.fastq reads

 1. Trimming and Adaptor Sequencing
  * Remove any Illumina Sequencing Adaptors, poly-A tails, sequence quality score <20, etc
 2. Overlapping of Read Libraries
  * We know our reads are 250bp PE with sequence overlaps, you might want to remove this step if yours are not.
 3. Assembly of Reads
  * We find that SPAdes gives good results, you may also like to try IDBA and Velvet-SC both of which can assemble SC.
 4. Assembly Statistics
  * Nice assembly statistics
 5. Read Mapping
  * This can help with coverage information etc, but we will be using it mostly for "blobology" (see step 7)
 6. BLAST Report
  * Top hits to NCBI's 'nt' database using 'megablast'
 7. BLOBTOOLS
  * GC/Coverage plots with taxonomy information to look for contamination.
 8. ?

# Downstream Analyses

 1. Contig Integrator for Sequence Assembly - [CISA](http://sb.nhri.org.tw/CISA/en/CISA)
 2. Protocol for fully automated Decontamination of Genomes - [ProDeGe](http://www.nature.com/ismej/journal/v10/n1/full/ismej2015100a.html)
 3. "Completeness" - Not convinced on these for SAGs
  * CEGMA?
  * BUSCO?

# Workflow Programs Required
 * [pigz](http://zlib.net/pigz/) - Parallel GZIP

 1. [Trim Galore!](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
  1. [cutadapt](https://cutadapt.readthedocs.org/en/stable/)
  2. [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
 2. [PEAR](http://sco.h-its.org/exelixis/web/software/pear/doc.html)
 3. [SPAdes](http://bioinf.spbau.ru/en/spades)
 4. [QUAST](http://bioinf.spbau.ru/quast)
 5. Mapping
  1. [BWA](https://github.com/lh3/bwa)
  2. [Samtools](http://www.htslib.org/)
 6. BLAST
  1. [BLAST+ executables](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) - megablast
  2. [NCBI 'nt' database](ftp://ftp.ncbi.nlm.nih.gov/blast/db/) - nt.*.tar.gz
  3. [NCBI Taxonomy dump](ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/) - taxdump.tar.gz
 7. [blobtools](https://github.com/DRL/blobtools)
 8. ?



# Dependencies
 1. [pigz](http://zlib.net/pigz/) - Parallel GZIP
 2. tee - GNU Core
 3. time - *nix Core
