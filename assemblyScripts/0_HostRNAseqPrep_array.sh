#!/bin/bash

# this script takes the raw transcriptome data from a separate experiment, trims illumina adapters,
	# and filters out endosymbiont contamination based on available Blattabacterium genomes
	# the bbsplit build contains blatta genomes, so I am keeping anything that does NOT align

ml BBMap/39.01-GCC-11.3.0

file=$(awk "NR==${SLURM_ARRAY_TASK_ID}" RNASamples.txt)

################################################################################################
##### Trim Adapters
################################################################################################
bbduk.sh in=RawData/"$file"_R1.fastq.gz \
out=RNAseqs/noadapt"$file"_R1.fastq.gz \
ref=/home/rld91076/adapters.fa \
ktrim=r k=23 mink=11 hdist=1 tpe tbo overwrite=true interleaved=f

bbduk.sh in=RawData/"$file"_R2.fastq.gz \
out=RNAseqs/noadapt"$file"_R2.fastq.gz \
ref=/home/rld91076/adapters.fa \
ktrim=r k=23 mink=11 hdist=1 tpe tbo overwrite=true interleaved=f

################################################################################################
##### Repair files
################################################################################################
repair.sh repair=t overwrite=true \
in1=RNAseqs/noadapt"$file"_R1.fastq.gz \
in2=RNAseqs/noadapt"$file"_R2.fastq.gz \
out=RNAseqs/rep"$file"_R1.fastq.gz \
out2=RNAseqs/rep"$file"_R2.fastq.gz \
outs=Unpaired/unpaired_"$file"_reads.fastq.gz

################################################################################################
##### Remove Blatta
################################################################################################
bbsplit.sh build=2 -Xmx100g requirecorrectstrand=f \
ambiguous=all ambiguous2=all qtrim=lr maxindel=20 ssao=t machineout=t \
in1=RNAseqs/rep"$file"_R1.fastq.gz \
in2=RNAseqs/rep"$file"_R2.fastq.gz \
outu1=BBSplit/Unmapped/noblatta"$file"_R1.fastq.gz \
outu2=BBSplit/Unmapped/noblatta"$file"_R2.fastq.gz \
2>&1 | tee BBSplit/results_"$file".txt

################################################################################################
##### Repair files again
################################################################################################
repair.sh repair=t overwrite=true \
in1=BBSplit/Unmapped/noblatta"$file"_R1.fastq.gz \
in2=BBSplit/Unmapped/noblatta"$file"_R2.fastq.gz \
out=PamRNAseqs/rep2"$file"_R1.fastq.gz \
out2=PamRNAseqs/rep2"$file"_R2.fastq.gz \
outs=Unpaired/stillunpaired_"$file"_reads.fastq.gz