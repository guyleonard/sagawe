# Single Cell Assembly Workflow

This is a suggested workflow for the assembly of MDA Single Cell MiSeq Illumina libraries. It is contained within a bash script that will execute all the programs sequentially.
The directory structure is assumed (read: adapted to our local situation) to be:

* WD/Sample_1/raw_illumina_reads/*.fastq
* ...
* WD/Sample_n/raw_illumina_reads/*.fastq

The script will iterate over each folder one by one...

## SAG Assembly Workflow in Steps

This description is mostly program agnostic, however, we have included below (and in the script) the programs we have found to be useful during our analyses and projects. You are welcome to adapt the workflow to inlude your favourite<sup>TM</sup> programs. 

1. Overlapping of Read Libraries<sup>[1](#footnote1)</sup>
  * We know our reads are 250bp PE with sequence overlaps, you might want to remove/adapt this step if yours are not.
2. Trimming and Adaptor Sequencing
  * Remove any Illumina Sequencing Adaptors, poly-A tails, sequence quality score <20, etc
3. Assembly of Reads
  * We find that SPAdes gives good results, you may also like to try [IDBA-UD](http://i.cs.hku.hk/~alse/hkubrg/projects/idba_ud/index.html) and [Velvet-SC](http://bix.ucsd.edu/projects/singlecell/) both of which can assemble SCs.
4. Assembly Statistics
  * Nice assembly statistics, e.g. N50, contig/scaffold lengths, etc.
5. Read Mapping
  * This can help with coverage information etc, but we will be using it mostly for "blobology" (see step 7)
6. BLAST Report
  * Top hits to NCBI's 'nt' database using 'megablast', beware false hits as with any BLAST search.
7. BLOBTOOLS
  * GC/Coverage plots with taxonomy information to look for contamination.
8. Genome 'Completeness' Test
9. Genome 'Completeness' Test
10. MultiQC - Aggregate results from bioinformatics analyses across many samples into a single report
11. ?

#### Other Dependencies
1. [pigz](http://zlib.net/pigz/) - Parallel GZIP
2. tee - GNU Core
3. time - *nix Core

## Assembly Downstream/Other Analyses

1. ESOM?
2. Contig Integrator for Sequence Assembly - [CISA](http://sb.nhri.org.tw/CISA/en/CISA)
3. Protocol for fully automated Decontamination of Genomes - [ProDeGe](http://www.nature.com/ismej/journal/v10/n1/full/ismej2015100a.html)
4. [CheckM](https://ecogenomics.github.io/CheckM/) - CheckM provides a set of tools for assessing the quality of genomes recovered from isolates, single cells, or metagenomes.
  * Looks like a nice tool, huge selection of options though and bac/arch oriented but can do ~euk.   
5. 

## Gene Prediction Workflow

I have used ansible to install the dependencies for this workflow. I have also included a method to download the repbase libraries using my password - however it is encrypted within the playbook, so it won't work for you. You will have to create your own ansible vault with this format

Also, rmblast won't currently download with Ansible 2.1.1.0 as there's something up with ftp downloads, so you will have do download it yourself to the .source dir.!?
    ---
    repbase_password: password

You can call the playbook to install like this:

    ansible-playbook install_gene_prediction_dependencies.yaml --sudo -K -c local -i "localhost," --ask-vault-pass

There are also tags so you can install one or many components in a go:

    ansible-playbook install_gene_prediction_dependencies.yaml --sudo -K -c local -i "localhost," --ask-vault-pass --tags repbase,hmmer

# Footnotes
<a name="footnote1">1</a>: Originally I had steps 1 and 2 in the reverse order (trimming and then overlapping), on some read libraries this caused issues (I think where paired reads would become unordered and so overlapping would not run), however I don't think this is the problem with the order. Reads should be overlapped first, prior to trimming, as we should end up with a set of longer reads - due to the better quality scores over-riding the lower qualities within the overlapped areas, where this would not have happened if the lower quality reads were already trimmed.
 
