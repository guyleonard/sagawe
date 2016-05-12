#!/bin/bash
set -x
set -e

start_dir=$(pwd)

# Setup environment variables
update_path () {
  new_dir=$1
  export PATH=${PATH:-$new_dir}
  if [[ ! "$PATH" =~ (^|:)"${new_dir}"(:|$) ]]; then
    export PATH=${new_dir}:${PATH}
  fi 
}

# Update Dependencies
sudo apt-get update -q
sudo apt-get install -y -q build-essential autoconf automake libtool python-setuptools python-dev python-pip pigz unzip default-jdk default-jre libfreetype6-dev libpng-dev pkg-config

# Make an install location
if [ ! -d 'build' ]; then
  mkdir build
fi
cd build
build_dir=$(pwd)

# Build all the things
cd $build_dir


# PEAR
echo "Downloading PEAR"
git clone https://github.com/xflouris/PEAR.git
cd PEAR
echo "Building PEAR"
./autogen.sh
./configure
make
sudo make install
echo "Done PEAR"
cd $build_dir
#


# Trim Galore!
echo "Trim Galore: cutadapt"
## cutadapt
sudo pip install --upgrade cutadapt

## FastQC
echo "Trim Galore: Downloading FastQC"
wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip
unzip fastqc_v0.11.5.zip
cd FastQC
fastqc_dir=$(pwd)
chmod 755 fastqc
echo "Trim Galore: Linking FastQC"
sudo ln -f -s $fastqc_dir/fastqc /usr/local/bin/fastqc
cd $build_dir

## Trim Galore
echo "Trim Galore: Downloading Trim Galore"
wget http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/trim_galore_v0.4.1.zip
unzip trim_galore_v0.4.1.zip
cd trim_galore_zip
trim_galore_dir=$(pwd)
chmod 755  trim_galore
echo "Trim Galore: Linking FastQC"
sudo ln -f -s $trim_galore_dir/trim_galore /usr/local/bin/trim_galore
cd $build_dir


# SPAdes
echo "SPAdes: Downloading SPAdes"
wget http://spades.bioinf.spbau.ru/release3.7.1/SPAdes-3.7.1-Linux.tar.gz
tar zxvf SPAdes-3.7.1-Linux.tar.gz
cd SPAdes-3.7.1-Linux/bin/
spades_dir=$(pwd)
echo "SPAdes: Updating PATH with SPAdes binary location"
update_path ${spades_dir}
cd $build_dir


# QUAST
echo "QUAST: Downloading"
wget https://sourceforge.net/projects/quast/files/quast-4.0.tar.gz
tar zxvf quast-4.0.tar.gz
cd quast-4.0
quast_galore_dir=$(pwd)
echo "QUAST: installing matplotlib"
sudo pip install --upgrade matplotlib
echo "QUAST: Testing to pre-compile"
python quast.py --test
update_path ${quast_dir}
cd $build_dir


# bwa
echo "BWA: Cloning from git"
git clone https://github.com/lh3/bwa.git
cd bwa
bwa_dir=$(pwd)
echo "BWA: Make"
make
echo "BWA: Installing"
sudo ln -s -f $bwa_dir/bwa /usr/local/bin/bwa_galore
cd $build_dir


# samtools
echo "samtools: Downloading"
wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2
bunzip2 samtools-1.3.1.tar.bz2
tar xvf samtools-1.3.1.tar
cd samtools-1.3.1
samtools_dir=$(pwd)
echo "samtools: Installing"
./configure
make
sudo make install
cd $build_dir


# blast
## blast+ executables
sudo "blast: blast+"
sudo apt-get install ncbi-blast+

## 'nt' database
#echo "blast: downloading 'nt' database"
#cd
#mkdir blast
#mkdir blast/nt
#cd blast/nt
#wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.*
#gunzip *.gz

## taxonomy db
echo "blast: downloading taxdump db"
cd
mkdir blast/taxonomy
cd blast/taxonomy
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar zxvf taxdump.tar.gz
cd $build_dir


## blobtools
echo "blobtools: Cloning github repository"
git clone https://github.com/DRL/blobtools.git
cd blobtools
blobtools_dir=$(pwd)
echo "blobtools: installing matplotlib"
#sudo pip install --upgrade matplotlib
# already required by QUAST
echo "blobtools: installing docopt"
sudo pip install docopt
update_bath blobtools_dir
cd $build_dir


# CEGMA
## geneid
echo "CEGMA: downloading geneid"
wget ftp://genome.crg.es/pub/software/geneid/geneid_v1.4.4.Jan_13_2011.tar.gz
cd geneid
geneid_dir=$(pwd)
sudo ln -f -s $geneid_dir/geneid /usr/local/bin/geneid

## genewise
echo "CEGMA: installing genewise"
sudo apt-get install wise

## hmmer
cd $build_dir
echo "CEGMA: installing hmmer"
wget http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz
tar zxvf hmmer-3.1b2-linux-intel-x86_64.tar.gz
cd hmmer-3.1b2-linux-intel-x86_64
hmmer_dir=$(pwd)
./configure
make
sudo make install
cd $build_dir


## CEGMA
echo "CEGMA: Downloading CEGMA v2.5
wget http://korflab.ucdavis.edu/datasets/cegma/CEGMA_v2.5.tar.gz
cd CEGMA_v2.5
cegma_dir=$(pwd)
make
cd $build_dir


#BUSCO
## Augustus
echo "BUSCO: downloading Augustus"
wget http://bioinf.uni-greifswald.de/augustus/binaries/augustus.current.tar.gz
tar zxvf augustus.current.tar.gz
cd augustus-3.2.2
augustus_dir=$(pwd)
echo "BUSCO: AUGUSTUS: Boost C++ & zlib & bamtools"
sudo apt-get install libboost-iostreams-dev libboost-graph-dev zlib1g-dev libgsl0-dev bamtools libbamtools-dev
echo "BUSCO: AUGUSTUS: Installing"
make
sudo make install
cd $build_dir

## Emboss Tools
echo "BUSCO: Installing EMBOSS Tools"
sudo apt-get install emboss

## BUSCO
echo "BUSCO: Downloading BUSCO"
wget http://busco.ezlab.org/files/BUSCO_v1.2.tar.gz
tar zxvf BUSCO_v1.2.tar.gz
cd BUSCO_v1.2
busco_dir=$(pwd)
cd $build_dir

## BUSCO DB
echo "BUSCO: Downloading Eukaryota DB"
cd
mkdir busco && cd busco
busco_db_dir=$(pwd)
wget http://busco.ezlab.org/files/eukaryota_buscos.tar.gz
tar zxvf eukaryota_buscos.tar.gz
cd $build_dir


# MultiQC
echo "MultiQC: Installing"
sudo pip install multiqc


# sometimes ld libraries are not linked, update them now
sudo ldconfig

echo "Done!"

set +x
set +e
