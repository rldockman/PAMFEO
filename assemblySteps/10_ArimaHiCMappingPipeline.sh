#!/bin/bash

##########################################################################################
# Notes
##########################################################################################
# followed this pipeline: https://github.com/ArimaGenomics/mapping_pipeline

##########################################################################################
# Modules
##########################################################################################
module load SAMtools/1.16.1-GCC-11.3.0
module load picard/2.27.5-Java-15

##########################################################################################
# Script
##########################################################################################
perl otherScripts/two_read_bam_combiner.pl filt_f.bam filt_r.bam samtools 10 | samtools view -@ 32 -bS -t filteredcontigs.fa.fai - | samtools sort -@ 32 -o combined.bam -

java -Xmx32G -Djava.io.tmpdir=temp/ -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups INPUT=combined.bam OUTPUT=rg.bam \
ID=primary LB=primary SM=PAMFEO PL=ILLUMINA PU=none

java -Xmx32G -XX:-UseGCOverheadLimit -Djava.io.tmpdir=temp/ -jar $EBROOTPICARD/picard.jar MarkDuplicates \
INPUT=rg.bam OUTPUT=PAMFEO_primary.bam METRICS_FILE=picardmetrics_primary.txt TMP_DIR=tempdir \
ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE

samtools index -@ 32 PAMFEO_primary.bam

perl otherScripts/get_stats.pl PAMFEO_primary.bam > stats/PAMFEO_primary.picardBam.stats