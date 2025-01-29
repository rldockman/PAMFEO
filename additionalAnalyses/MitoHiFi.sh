#!/bin/bash

# First create singularity image from docker
apptainer pull docker://ghcr.io/marcelauliano/mitohifi:master

# If you need to download a mitochondrial reference:
singularity exec ~/Genome/mitohifi_master.sif \
findMitoReference.py --species "Periplaneta americana" --outfolder ~/Genome/mitohifi --min_length 14000

# Used original hifiasm output fasta file, since I found no matches in filtered genome (all were removed)
singularity exec ~/Genome/mitohifi_master.sif \
mitohifi.py -c ~/Genome/hifiasm/pamerasm_p.fa -f NC_016956.1.fasta -g NC_016956.1.gb -t 16 -o 5