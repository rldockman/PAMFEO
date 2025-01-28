#!/bin/bash

##########################################################################################
# Notes
##########################################################################################
# ran separate script per tissue: HgT, MgT, FB

##########################################################################################
# Modules
##########################################################################################
module load SAMtools/1.16.1-GCC-11.3.0

##########################################################################################
# Script
##########################################################################################
### fat body
samtools merge -@ 16 -o merged_FB.bam \
filt_maskedXT1_FB1_sort.bam filt_maskedXT1_FB2_sort.bam filt_maskedXT1_FB3_sort.bam filt_maskedXT1_FB5_sort.bam \
filt_maskedXT5_FB1_sort.bam filt_maskedXT5_FB2_sort.bam filt_maskedXT5_FB3_sort.bam filt_maskedXT5_FB6_sort.bam

samtools coverage -o stats/FBcov.tsv merged_FB.bam

### hindgut tissue
samtools merge -@ 16 -o merged_HgT.bam \
filt_maskedXT1_HgT1_sort.bam filt_maskedXT1_HgT2_sort.bam filt_maskedXT1_HgT3_sort.bam filt_maskedXT1_HgT5_sort.bam \
filt_maskedXT5_HgT1_sort.bam filt_maskedXT5_HgT2_sort.bam filt_maskedXT5_HgT3_sort.bam filt_maskedXT5_HgT6_sort.bam

samtools coverage -o stats/HgTcov.tsv merged_HgT.bam

### midgut tissue
samtools merge -@ 16 -o merged_MgT.bam \
filt_maskedXT1_MgT1_sort.bam filt_maskedXT1_MgT2_sort.bam filt_maskedXT1_MgT3_sort.bam filt_maskedXT1_MgT5_sort.bam \
filt_maskedXT5_MgT1_sort.bam filt_maskedXT5_MgT2_sort.bam filt_maskedXT5_MgT3_sort.bam filt_maskedXT5_MgT6_sort.bam

samtools coverage -o stats/MgTcov.tsv merged_MgT.bam