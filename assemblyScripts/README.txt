These are the scripts used to generate the scaffolded primary genome assembly from HiFi Pacbio long reads and HiC short reads.
Each file contains multiple scripts, divided by hashes, as they were submitted to UGA's high performance computing cluster.
    I prefered to check the output in between steps and needed to do so for some manual filtering steps, but with enough resources allocated, some steps can be combined.

0_HostRNAseqPrep_array.sh
    uses JGI's BBmap/bbsplit/bbduk programs to trim and sort paired-end raw illumina RNAseq reads that are used in the filtering process

1_PAMFEO_AssemblyandFiltering_collapsed.sh
    first set of scripts that: 
    generate initial assembly with hifiasm,
    calculate Pacbio coverage statistics, 
    filter out contamination via kraken2, 
    align and assess RNAseq coverage of each contig,
    filter contigs from assembly to consensus set for scaffolding

2_PAMFEO_HiCAlignmentandScaffolding_collapsed.sh
    second set of scripts that:
    follow the Arima Genomics pipeline for HiC read alignment,
    scaffold the contig assembly and chromatin contact reads via YaHS,
    finalize the assembly after curation in Juicebox and telomere identification (TIDK.sh),
    calculate genome BUSCO score on final assembly

Tools used:
BBMap/39.01
BamTools/2.5.1
BLAST+/2.12.0
pigz/2.7
hifiasm/0.19.6
SeqKit/2.5.1
Kraken2/2.1.3
RepeatMasker/4.1.5
HISAT2/3n-20201216
SAMtools/1.16.1
BWA/0.7.17
picard/2.27.5
BUSCO/5.5.0
BEDTools/2.30.0
YaHS/1.1
Juicebox/2.20.00
