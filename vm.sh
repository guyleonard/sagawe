sudo apt update
sudo apt upgrade
sudo apt dist-upgrade

sudo apt install python-pip unzip openjdk-11-jdk pigz zlib1g-dev libncurses5-dev libbz2-dev liblzma-dev r-base-core cpanminus

sudo -H pip install cutadapt

mkdir packages && cd packages

wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip
unzip && fastqc_v0.11.8.zip
chmod +x fastqc
PATH=$PATH:/home/ubuntu/packages/FastQC

curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.4.tar.gz -o trim_galore.tar.gz
tar xvzf trim_galore.tar.gz
PATH=$PATH:/home/ubuntu/packages/TrimGalore-0.6.4

wget https://sourceforge.net/projects/bbmap/files/BBMap_38.67.tar.gz
tar zxvf BBMap_38.67.tar.gz
PATH=$PATH:/home/ubuntu/packages/bbmap

wget http://cab.spbu.ru/files/release3.13.1/SPAdes-3.13.1-Linux.tar.gz
tar -xzf SPAdes-3.13.1-Linux.tar.gz
PATH=$PATH:/home/ubuntu/packages/SPAdes-3.13.1-Linux

wget https://downloads.sourceforge.net/project/quast/quast-5.0.2.tar.gz
tar -xzf quast-5.0.2.tar.gz
PATH=$PATH:/home/ubuntu/packages/quast-5.0.2

wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.9.0+-x64-linux.tar.gz
tar zxvf ncbi-blast-2.9.0+-x64-linux.tar.gz
PATH=$PATH:/home/ubuntu/packages/ncbi-blast-2.9.0+/bin

wget https://github.com/DRL/blobtools/archive/blobtools_v1.1.1.tar.gz
tar zxvf blobtools_v1.1.1.tar.gz
PATH=$PATH:/home/ubuntu/packages/blobtools-blobtools_v1.1.1
XXXXXXX

git clone https://github.com/lh3/bwa.git
cd bwa; make
PATH=$PATH:/home/ubuntu/packages/bwa

wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
bunzip2 samtools-1.9.tar.bz2 && tar xvf samtools-1.9.tar
cd samtools-1.9 && ./configure
make && sudo make install

git clone https://gitlab.com/ezlab/busco.git
sudo python setup.py install
PATH=$PATH:/home/ubuntu/packages/busco/scripts
mkdir databases && cd databases
wget http://busco.ezlab.org/v2/datasets/bacteria_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/eukaryota_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/fungi_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/protists_ensembl.tar.gz
for i in *.gz; do tar zxvf $i; done

