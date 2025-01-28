#!/bin/bash

##########################################################################################
# Notes
##########################################################################################
# Make array with RNA samples so each can be done separately

##########################################################################################
# Modules
##########################################################################################
module load HISAT2/3n-20201216-gompi-2022a
module load SAMtools/1.16.1-GCC-11.3.0

##########################################################################################
# Script
##########################################################################################
file=$(awk "NR==${SLURM_ARRAY_TASK_ID}" RNASamples.txt)

hisat2 -p 32 -x hisat2/maskedgenome -k 3 --dta --no-unal \
-1 Transcriptome/PamRNAseqs/rep2"$file"_R1.fastq.gz -2 Transcriptome/PamRNAseqs/rep2"$file"_R2.fastq.gz \
--rg-id "$file" --rg SM:"$file" --summary-file hisat2/hisat2summary_"$file".txt | \
samtools view -@32 -Sb - | \
samtools sort -@32 --write-index - -o filt_masked"$file"_sort.bam##idx##filt_masked"$file"_sort.bam.bai