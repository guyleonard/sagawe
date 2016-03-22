# Single Cell Assembly Workflow

This is a suggested workflow for the assembly of MDA Single Cell MiSeq Illumina libraries.

## Workflow in Steps

 1. Trimming and Adaptor Sequencing
 2. Overlapping of Read Libraries
 3. Assembly of Reads
 4. Assembly Statistics
 5. Read Mapping
 6. BLAST Report
 7. BLOBTOOLS
 8. ?

## Workflow Programs
 1. [Trim Galore!](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
 ..* [cutadapt](https://cutadapt.readthedocs.org/en/stable/)
 ..* [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
 2. [PEAR](http://sco.h-its.org/exelixis/web/software/pear/doc.html)
 3. [SPAdes](http://bioinf.spbau.ru/en/spades)
 4. [QUAST](http://bioinf.spbau.ru/quast)
 5. Mapping
 ..* [BWA](https://github.com/lh3/bwa)
 ..* [Samtools](http://www.htslib.org/)
 6. BLAST
 ..* [BLAST+ executables](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) - megablast
 ..* [NCBI 'nt' database](ftp://ftp.ncbi.nlm.nih.gov/blast/db/) - nt.*.tar.gz
 ..* [NCBI Taxonomy dump](ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/) - taxdump.tar.gz
 7. [blobtools](https://github.com/DRL/blobtools)
 8. ?
