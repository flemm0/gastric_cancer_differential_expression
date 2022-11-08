#!/bin/bash

##################################################################

# make directories

cd /scratch/
mkdir FASTQS REFS Tophat-BAMS Cuffdiff-STATS Cuffdiff-noSTATS


##################################################################


# download SRA sequences as fastq files

cd FASTQS

for (( i = 70; i <= 77; i++ ))
do
	nohup fasterq-dump -S SRR5855$i &
done


##################################################################

# check read counts

for (( i = 70; i <= 77; i++ ))
do
	echo "Reads for SRR5855${i}_1:"
	wc -l SRR5855${i}_1.fastq | awk '{x=$1/4; print x}'
	echo "Reads for SRR5855${i}_1:"
	wc -l SRR5855${i}_2.fastq | awk '{x=$1/4; print x}'
done


##################################################################

# check that read 1 and read 2 are not the same

for (( i = 70; i <= 77; i++ ))
do
	echo "SRR5855${i}_1:"
	head -n2 SRR5855${i}_1.fastq
	echo "SRR5855${i}_2:"
	head -n2 SRR5855${i}_1.fastq
	printf '\n'
done


##################################################################

# check cpu architecture

# lscpu

### server cpu architecture is x86_64


##################################################################

# download human reference genome fasta and gtf files from ensembl

cd ../REFS

nohup wget https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz &

nohup wget https://ftp.ensembl.org/pub/current_gtf/homo_sapiens/Homo_sapiens.GRCh38.108.chr.gtf.gz &

gunzip *.gz


##################################################################

# build index

nohup bowtie2-build --threads 12 Homo_sapiens.GRCh38.dna.primary_assembly.fa Homo_sapiens.GRCh38.dna.primary_assembly > bowtie2.nohup.out &


##################################################################

# preliminary alignment to whole genome and whole transcriptome
## will be aligning SRR85571 (gastric normal tissue) and SRR585573 (gastric cancer tissue cell culture)


cd ../FASTQS

# align SRR585571 to genome

nohup tophat2 \
	-p 8 \
	-G ../REFS/Homo_sapiens.GRCh38.108.chr.gtf \
	-o ./sample1_5571_genomeAlign \
	../REFS/Homo_sapiens.GRCh38.dna.primary_assembly \
	SRR585571_1.fastq,SRR585571_2.fastq > \
	sample1_5571.genome.nohup.out &


# align SRR585571 to transcriptome

nohup tophat2 \
	-p 8 \
	-G ../REFS/Homo_sapiens.GRCh38.108.chr.gtf \
	-T \
	-o ./sample1_5571_transcriptomeAlign \
	../REFS/Homo_sapiens.GRCh38.dna.primary_assembly \
	SRR585571_1.fastq,SRR585571_2.fastq > \
	sample1_5571.transcriptome.nohup.out &


# align SRR585573 to genome

nohup tophat2 \
	-p 12 \
	-G ../REFS/Homo_sapiens.GRCh38.108.chr.gtf \
	-o ./sample2_5573_genomeAlign \
	../REFS/Homo_sapiens.GRCh38.dna.primary_assembly \
	SRR585573_1.fastq,SRR585573_2.fastq > \
	sample2_5573.genome.nohup.out &


# align SRR585573 to transcriptome

nohup tophat2 \
	-p 12 \
	-G ../REFS/Homo_sapiens.GRCh38.108.chr.gtf \
	-T \
	-o ./sample2_5573_transcriptomeAlign \
	../REFS/Homo_sapiens.GRCh38.dna.primary_assembly \
	SRR585573_1.fastq,SRR585573_2.fastq > \
	sample2_5573.transcriptome.nohup.out &


##################################################################

# align remaining reads to whole genome

for i in 70 72 74 75 76 77
do
	nohup tophat2 \
		-p 12 \
		-G ../REFS/Homo_sapiens.GRCh38.108.chr.gtf \
		-o ./5855${i}_genomeAlign \
		../REFS/Homo_sapiens.GRCh38.dna.primary_assembly \
		SRR5855${i}_1.fastq,SRR5855${i}_2.fastq > \
		5855${i}.genome.nohup.out &
done


##################################################################

# create symlinks for bamfiles in FASTQS directory to access within Tophat-BAMS directory

cd ../Tophat-BAMS

ln -s /scratch/wuflemmi-F22/FASTQS/585570_genomeAlign/accepted_hits.bam normal1
ln -s /scratch/wuflemmi-F22/FASTQS/sample1_585571_genomeAlign/accepted_hits.bam normal1
ln -s /scratch/wuflemmi-F22/FASTQS/585572_genomeAlign/accepted_hits.bam cancer_tissue1
ln -s /scratch/wuflemmi-F22/FASTQS/sample2_585573_genomeAlign/accepted_hits.bam cancer_tissue2
ln -s /scratch/wuflemmi-F22/FASTQS/585574_genomeAlign/accepted_hits.bam cancer_tissue3 
ln -s /scratch/wuflemmi-F22/FASTQS/585575_genomeAlign/accepted_hits.bam cancer_cl1
ln -s /scratch/wuflemmi-F22/FASTQS/585576_genomeAlign/accepted_hits.bam cancer_cl2
ln -s /scratch/wuflemmi-F22/FASTQS/585577_genomeAlign/accepted_hits.bam cancer_cl3

