sudo apt update
sudo apt upgrade
sudo apt dist-upgrade

sudo apt install python-pip unzip openjdk-11-jdk pigz zlib1g-dev libncurses5-dev libbz2-dev liblzma-dev r-base-core cpanminus libboost-iostreams-dev zlib1g-dev libbamtools-dev libssl-dev libcurl3-dev cmake libjsoncpp-dev python3-dev python3-pip libxml2-dev wise

sudo -H pip install cutadapt

mkdir packages && cd packages

wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip
unzip && fastqc_v0.11.8.zip
cd FastQC
chmod +x fastqc
PATH=$PATH:/home/ubuntu/packages/FastQC
cd ..

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
sudo cpanm JSON
mkdir databases && cd databases
perl ~/packages/ncbi-blast-2.9.0+/bin/update_blastdb.pl taxdb
perl ~/packages/ncbi-blast-2.9.0+/bin/update_blastdb.pl nt

wget https://github.com/DRL/blobtools/archive/blobtools_v1.1.1.tar.gz
tar zxvf blobtools_v1.1.1.tar.gz
sudo pip3 install --upgrade pip
sudo -H pip3 install docopt
sudo python setup.py install
PATH=$PATH:/home/ubuntu/packages/blobtools-blobtools_v1.1.1
cd ~/databases
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar zxf taxdump.tar.gz  nodes.dmp names.dmp
blobtools nodesdb --nodes ./nodes.dmp --names ./names.dmp

git clone https://github.com/lh3/bwa.git
cd bwa; make
PATH=$PATH:/home/ubuntu/packages/bwa
cd ..

git clone https://github.com/samtools/htslib.git
cd htslib
autoheader && autoconf
./configure && make
sudo make install
cd ..

git clone https://github.com/samtools/bcftools.git
cd bcftools
autoheader &&autoconf
./configure && make
sudo make install
cd 

export TOOLDIR=/home/ubuntu/packages

git clone https://github.com/Gaius-Augustus/Augustus.git
cd Augustus
make && sudo make install
export AUGUSTUS_CONFIG_PATH=/home/ubuntu/packages/Augustus/config
cd ..

git clone https://gitlab.com/ezlab/busco.git
cd busco
sudo python setup.py install
PATH=$PATH:/home/ubuntu/packages/busco/scripts
cd ../databases
wget http://busco.ezlab.org/v2/datasets/bacteria_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/eukaryota_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/fungi_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/protists_ensembl.tar.gz
for i in *.gz; do tar zxvf $i; done
cd ..

git clone git://github.com/pezmaster31/bamtools.git
cd bamtools && mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=/home/ubuntu/packages/bamtools ..
make && sudo make install
cd ..

wget http://eddylab.org/software/hmmer/hmmer.tar.gz
tar zxvf hmmer.tar.gz
cd hmmer-3.2.1
./configure && make
cd ..

git clone https://github.com/TGAC/KAT.git
cd KAT
./build_boost.sh
./autogen.sh
sudo -H pip3 install numpy scipy matplotlib sphinx tabulate
./configure
make && sudo make install
cd ..

git clone https://github.com/schatzlab/genomescope.git
PATH=$PATH:/home/ubuntu/packages/genomescope

git clone https://github.com/KamilSJaron/smudgeplot
cd smudgeplot
sudo Rscript -e 'install.packages("devtools")'
sudo Rscript -e 'install.packages("argparse")'
sudo Rscript -e 'install.packages("viridis")'
sudo Rscript install.R
sudo install -C exec/smudgeplot.py /usr/local/bin
sudo install -C exec/smudgeplot_plot.R /usr/local/bin
cd ..

wget https://github.com/gmarcais/Jellyfish/releases/download/v2.3.0/jellyfish-2.3.0.tar.gz
tar zxvf jellyfish-2.3.0.tar.gz
cd jellyfish-2.3.0
./configure && make

sudo -H pip install multiqc

wget https://github.com/KorfLab/CEGMA_v2/archive/v2.5.tar.gz
tar zxvf v2.5.tar.gz
make
PATH=$PATH:/home/ubuntu/packages/CEGMA_v2-2.5/bin
export CEGMA=/home/ubuntu/packages/CEGMA_v2-2.5
export PERL5LIB="$PERL5LIB:$CEGMA/lib"
export WISECONFIGDIR=/usr/share/wise/

wget ftp://genome.crg.es/pub/software/geneid/geneid_v1.4.4.Jan_13_2011.tar.gz
tar zxvf geneid_v1.4.4.Jan_13_2011.tar.gz
cd geneid
make
PATH=$PATH:/home/ubuntu/packages/geneid/bin
