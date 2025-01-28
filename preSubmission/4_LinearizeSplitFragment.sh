#!/bin/bash

##########################################################################################
# Notes
##########################################################################################
# linearizing genome sets each contig on one continuous line within multifasta
	# easy to grab desired contigs
# this script fragments contigs for contamination filtering

##########################################################################################
# Modules
##########################################################################################
module load SeqKit/2.5.1

##########################################################################################
# Script
##########################################################################################
cd PBfiltered
mkdir contig_ind # for individual contig .fa files
mkdir contig_frag # for fragmented contig .fa files

# Linearize genome so every contig sequence is one line long
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < PAMFEO_p.fa > int.fa
tail -n +2 int.fa > linear_PAMFEO.fa
rm int.fa

cd contig_ind

# Place every contig into its own file
cat ~/PBfiltered/linear_PAMFEO.fa| awk '{
        if (substr($0, 1, 1)==">") {filename=(substr($0,2) ".fa")}
        print $0 >> filename
        close(filename)
}'

# Make file of contig names for all .fa files in current directory
ls *.fa > files.txt 

cd ~/PBfiltered

# Use seqkit to split contigs into 100nt fragments that overlap by 50nt
for line in $(cat contig_ind/files.txt); do
seqkit sliding -s 50 -W 100 contig_ind/$line > contig_frag/frag"$line"
done

cd contig_frag
ls *.fa > fragfiles.txt