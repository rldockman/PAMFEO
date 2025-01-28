#!/bin/bash

##########################################################################################
# Notes
##########################################################################################
### I am starting partway through the process with the blocklist file generated already

##########################################################################################
# Modules
##########################################################################################
module load BamTools/2.5.1-GCC-10.2.0
module load BLAST+/2.12.0-gompi-2020b
module load pigz/2.7-GCCcore-10.2.0

##########################################################################################
# Script
##########################################################################################
for x in `ls *bam | sed 's/\.bam//g'`;
do
bamtools convert -format fastq -in ${x}.bam -out ${x}.fastq
done

for x in `ls *fastq | sed 's/\.fastq//g'`;
do
cat ${x}.fastq | paste - - - - | grep -v -f ${x}.blocklist -F | tr "\t" "\n" | pigz -p 40 --fast > ${x}.fcsfilt.fastq.gz
f=`cat ${x}.blocklist | wc -l` #number of adapter contaminated 
r1=`cat ${x}.fastq | wc -l` 
r2=`awk -v r1=$r1 'BEGIN{ans=r1/4; print ans}'` #number of ccs reads
p1=`awk -v n1=$r2 -v n2=$f 'BEGIN{ans=n2/n1*100; print ans}'` #proportion of adapter contaminated reads
r3=`awk -v r2=$r2 -v f=$f 'BEGIN{ans=r2-f; print ans}'` #number of reads retained
p2=`awk -v p1=$p1 'BEGIN{ans=100-p1; print ans}'` #proportion of reads retained
echo "For the" ${x} "dataset:" >>${x}.stats 
echo "" >>${x}.stats
echo "Number of ccs reads:" $r2 >>${x}.stats
echo "Number of adapter contaminated ccs reads:" $f "("$p1"% of total)" >>${x}.stats
echo "Number of ccs reads retained:" $r3 "("$p2"% of total)" >>${x}.stats
echo "" >>${x}.stats
echo "Finished on $(date)" >>${x}.stats
done