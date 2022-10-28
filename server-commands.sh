#!/bin/bash

# make directories

cd /scratch/
mkdir FASTQS REFS Tophat-BAMS Cuffdiff-STATS Cuffdiff-noSTATS


# download SRA sequences as fastq files

cd FASTQS

for (( i = 70; i <= 77; i++ ))
do
	nohup fasterq-dump -S SRR5855$i &
done

# check read counts

for (( i = 70; i <= 77; i++ ))
do
	echo "Reads for SRR5855${i}_1:"
	wc -l SRR5855${i}_1.fastq | awk '{x=$1/4; print x}'
	echo "Reads for SRR5855${i}_1:"
	wc -l SRR5855${i}_2.fastq | awk '{x=$1/4; print x}'
done


# check that read 1 and read 2 are not the same

for (( i = 70; i <= 77; i++ ))
do
	echo "SRR5855${i}_1:"
	head -n2 SRR5855${i}_1.fastq
	echo "SRR5855${i}_2:"
	head -n2 SRR5855${i}_1.fastq
	printf '\n'
done
