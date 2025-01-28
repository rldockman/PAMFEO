#!/bin/bash

##########################################################################################
# Notes
##########################################################################################
# I ran YaHS without splitting contigs since the hifiasm assembly was high quality to begin with
# running with default parameters made it messier without improving the quality

# Note: juicer_tools.2 did not work with this version of YaHS
	# had to download older .jar file, and point towards that directory to work properly

##########################################################################################
# Modules
##########################################################################################
module load YaHS/1.1-GCC-11.3.0
module load Juicebox/2.20.00
module load SAMtools/1.16.1-GCC-11.3.0

##########################################################################################
# Script
##########################################################################################
# limits scaffolding only to contigs > 100,000bp (param -l), do not split contigs (param --no-contig-ec)
yahs --no-contig-ec -l 100000 filteredcontigs.fa PAMFEO_primary.bed 

# Index scaffolds file and generate chromosome sizes
samtools faidx yahs.out_scaffolds_final.fa
cut -f 1,2 yahs.out_scaffolds_final.fa.fai > scaffolds_final.chrom.sizes

### To visualize:
# Obtain alignments with "juicer pre" from yahs with bin file, agp file, and the index from the original non-scaffolded genome fasta
(juicer pre yahs.out.bin yahs.out_scaffolds_final.agp filteredcontigs.fa.fai | sort -k2,2d -k6,6d -T ./ --parallel=32 -S100G | awk 'NF' > alignments_sorted.txt.part) && (mv alignments_sorted.txt.part alignments_sorted.txt)

# Generate .hic files
(java -jar -Xmx150G otherScripts/juicer_tools.1.9.9_jcuda.0.8.jar pre alignments_sorted.txt out.hic.part scaffolds_final.chrom.sizes) && (mv out.hic.part out.hic)

### For manual curation in Juicebox
# Run juicer pre command with -a flag
juicer pre -a -o out_JBAT yahs.out.bin yahs.out_scaffolds_final.agp filteredcontigs.fa.fai >out_JBAT.log 2>&1

# Use output to generate different .hic file before, that can be messed with
# Look at out_JBAT.log to find assembly length (halved)
java -Xmx150G -jar otherScripts/juicer_tools.1.9.9_jcuda.0.8.jar pre out_JBAT.txt out_JBAT.hic <(echo "assembly 1615942301")