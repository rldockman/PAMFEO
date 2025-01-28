#!/bin/bash

##########################################################################################
# Notes
##########################################################################################
# converted both haplotypes and primary assembly
# manually curated .tsv files afterwards based on PacBio coverage
	# I kept contigs with coverage between 6 and 30 (rd:i number)

##########################################################################################
# Script
##########################################################################################
# Convert to fasta file
awk '/^S/{print ">"$2;print $3}' hifiasm/PAMFEO.asm.hic.p_ctg.gfa > PAMFEO_p.fa
awk '/^S/{print ">"$2;print $3}' hifiasm/PAMFEO.asm.hic.hap1.p_ctg.gfa > hap/PAMFEO_hap1.fa
awk '/^S/{print ">"$2;print $3}' hifiasm/PAMFEO.asm.hic.hap2.p_ctg.gfa > hap/PAMFEO_hap2.fa

# Obtain hifiasm data for contig length and coverage
grep ^S hifiasm/PAMFEO.asm.hic.p_ctg.noseq.gfa > stats/PBcov_primary.tsv
grep ^S hifiasm/PAMFEO.asm.hic.hap1.p_ctg.noseq.gfa > stats/PBcov_hap1.tsv
grep ^S hifiasm/PAMFEO.asm.hic.hap2.p_ctg.noseq.gfa > stats/PBcov_hap2.tsv