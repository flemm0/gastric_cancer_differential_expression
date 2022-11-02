#Tools needed: sratoolkit, samtools, bedtools, gatk, tophat2, bowtie2
#Other files needed: human reference genome

########################
#install sratoolkit

cd ~/TRGN514/ngs-tools

wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.0/sratoolkit.3.0.0-centos_linux64.tar.gz

tar -zxvf sratoolkit.3.0.0-centos_linux64.tar.gz

rm sratoolkit.3.0.0-centos_linux64.tar.gz

echo 'export PATH=$PATH:/home/wuflemmi-F22/TRGN514/ngs-tools/sratoolkit.3.0.0-centos_linux64/bin' >> ~/.bashrc



########################
# install samtools

wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2

tar -jxvf samtools-1.16.1.tar.bz2

rm samtools-1.16.1.tar.bz2

cd samtools-1.16.1

./configure

make

echo 'export PATH=$PATH:/home/wuflemmi-F22/TRGN514/ngs-tools/samtools-1.16.1' >> ~/.bashrc

cd ..

#######################
# install bedtools2

wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools-2.30.0.tar.gz

tar -xzvf bedtools-2.30.0.tar.gz

make

rm bedtools-2.30.0.tar.gz

echo 'export PATH=$PATH:/home/wuflemmi-F22/TRGN514/ngs-tools/bedtools2/bin' >> ~/.bashrc



#######################
# install bowtie2

wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.5/bowtie2-2.4.5-linux-x86_64.zip/download

unzip bowtie2-2.4.5-linux-x86_64.zip

rm bowtie2-2.4.5-linux-x86_64.zip

echo 'export PATH=$PATH:/home/wuflemmi-F22/TRGN514/ngs-tools/bowtie2-2.4.5-linux-x86_64' >> ~/.bashrc



#######################
# install tophat2

wget http://ccb.jhu.edu/software/tophat/downloads/tophat-2.1.1.Linux_x86_64.tar.gz

tar -xzvf tophat-2.1.1.Linux_x86_64.tar.gz

rm tophat-2.1.1.Linux_x86_64.tar.gz

echo 'export PATH=$PATH:/home/wuflemmi-F22/TRGN514/ngs-tools/tophat-2.1.1.Linux_x86_64' >> ~/.bashrc



#######################
# install GATK

wget https://github.com/broadinstitute/gatk/releases/download/4.3.0.0/gatk-4.3.0.0.zip

unzip gatk-4.3.0.0.zip

rm gatk-4.3.0.0.zip

echo 'export GATK=/home/wuflemmi-F22/TRGN514/ngs-tools/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar' >> ~/.bashrc


source ~/.bashrc




#######################
# download human reference genome (GrCh38)

cd /scratch/wuflemmi-F22/REFS

wget ftp://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz

gunzip Homo_sapiens.GRCh38.dna.toplevel.fa.gz &
