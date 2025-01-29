# Sections indicate separate cluster submission scripts

# Modules used:
BWA/0.7.17
SAMtools/1.16.1
picard/2.27.5
BUSCO/5.5.0
BEDTools/2.30.0
YaHS/1.1
Juicebox/2.20.00

# Starting genome assembly
filteredcontigs.fa # primary assembly from hifiasm, filtered by PacBio coverage, RNA-seq coverage, and Kraken2 hits

# Finished genome assembly, submitted to NCBI
PAMFEO1_hifiasm_priV1.fa

##############################################################################
# 9a_IndexGenome.sh
##############################################################################
# index filtered assembly for BWA and samtools

ml BWA/0.7.17-GCCcore-11.3.0
ml SAMtools/1.16.1-GCC-11.3.0

bwa index filteredcontigs.fa -a bwtsw
samtools faidx filteredcontigs.fa

##############################################################################
# 9b_BWA_HiC_F.sh
##############################################################################
# align forward HiC reads to filtered assembly

ml BWA/0.7.17-GCCcore-11.3.0
ml SAMtools/1.16.1-GCC-11.3.0

bwa mem -t 64 -T 32 filteredcontigs.fa \
~/Genome/RawData/HiCPAMFEO_R1.fastq.gz|samtools view -b -@ 64 - > pam_f.bam 

samtools view -@ 64 -h pam_f.bam | perl ~/Genome/scripts/filter_five_end.pl | samtools view -@ 64 -Sb - > filt_f.bam

##############################################################################
# 9c_BWA_HiC_R.sh
##############################################################################
# align reverse HiC reads to filtered assembly

ml BWA/0.7.17-GCCcore-11.3.0
ml SAMtools/1.16.1-GCC-11.3.0

bwa mem -t 64 -T 32 filteredcontigs.fa \
~/Genome/RawData/HiCPAMFEO_R2.fastq.gz|samtools view -b -@ 64 - > pam_r.bam 

samtools view -@ 64 -h pam_r.bam | perl ~/Genome/scripts/filter_five_end.pl | samtools view -@ 64 -Sb - > filt_r.bam

##############################################################################
# 10_ArimaHiCMappingPipeline.sh
##############################################################################
# see https://github.com/ArimaGenomics/mapping_pipeline for details and perl scripts

ml SAMtools/1.16.1-GCC-11.3.0
ml picard/2.27.5-Java-15

cd $SLURM_SUBMIT_DIR

perl ~/Genome/scripts/two_read_bam_combiner.pl filt_f.bam filt_r.bam samtools 10 | samtools view -@ 32 -bS -t filteredcontigs.fa.fai - | samtools sort -@ 32 -o combined.bam -

java -Xmx32G -Djava.io.tmpdir=temp/ -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups INPUT=combined.bam OUTPUT=rg.bam \
ID=primary LB=primary SM=PAMFEO PL=ILLUMINA PU=none

java -Xmx32G -XX:-UseGCOverheadLimit -Djava.io.tmpdir=temp/ -jar $EBROOTPICARD/picard.jar MarkDuplicates \
INPUT=rg.bam OUTPUT=PAMFEO_primary.bam METRICS_FILE=picardmetrics_primary.txt TMP_DIR=tempdir \
ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE

samtools index -@ 32 PAMFEO_primary.bam

perl ~/Genome/scripts/get_stats.pl PAMFEO_primary.bam > PAMFEO_primary.bam.stats

##############################################################################
# 11_HiC_Coverage.sh
##############################################################################
# for data on HiC read coverage

module load SAMtools/1.16.1-GCC-11.3.0

samtools coverage -o HiCcov.tsv PAMFEO_primary.bam

##############################################################################
# 12_BUSCO_PreScaf.sh
##############################################################################
# busco

module load BUSCO/5.5.0-foss-2022a

busco -m genome -r -i filteredcontigs.fa -o buscoPreScaf -l insecta_odb10 --cpu 16 --download_path ~/Genome/busco_downloads

##############################################################################
# 13_BAM2BED.sh
##############################################################################
# convert .bam to .bed for scaffolding
 
ml BEDTools/2.30.0-GCC-11.3.0
ml SAMtools/1.16.1-GCC-11.3.0

bamToBed -i PAMFEO_primary.bam > PAMFEO_primary.bed
sort -k 4 PAMFEO_primary.bed > tmp && mv tmp PAMFEO_primary.bed

##############################################################################
# 14_YaHS_scaffolding.sh
##############################################################################
# scaffold genome based on HiC links with YaHS
# Note: juicer_tools.2 did not work for me. Had to download older .jar file, and point towards that directory to work properly

module load YaHS/1.1-GCC-11.3.0
module load Juicebox/2.20.00
module load SAMtools/1.16.1-GCC-11.3.0

# Run yahs with fasta file and BED file of HiC alignments
yahs --no-contig-ec -l 100000 filteredcontigs.fa PAMFEO_primary.bed 

# Index scaffolds file and generate chromosome sizes
samtools faidx yahs.out_scaffolds_final.fa
cut -f 1,2 yahs.out_scaffolds_final.fa.fai > scaffolds_final.chrom.sizes

### To visualize:
# Obtain alignments with "juicer pre" from yahs with bin file, agp file, and the index from the original non-scaffolded genome fasta
(juicer pre yahs.out.bin yahs.out_scaffolds_final.agp filteredcontigs.fa.fai | sort -k2,2d -k6,6d -T ./ --parallel=32 -S100G | awk 'NF' > alignments_sorted.txt.part) && (mv alignments_sorted.txt.part alignments_sorted.txt)

# Generate .hic files
(java -jar -Xmx150G ~/Genome/yahs/juicer_tools.1.9.9_jcuda.0.8.jar pre alignments_sorted.txt out.hic.part scaffolds_final.chrom.sizes) && (mv out.hic.part out.hic)

### For manual curation in Juicebox
# Run juicer pre command with -a flag
juicer pre -a -o out_JBAT yahs.out.bin yahs.out_scaffolds_final.agp filteredcontigs.fa.fai >out_JBAT.log 2>&1

# Use output to generate different .hic file before, that can be messed with
# Look at out_JBAT.log to find assembly length (halved)
java -Xmx150G -jar ~/Genome/yahs/juicer_tools.1.9.9_jcuda.0.8.jar pre out_JBAT.txt out_JBAT.hic <(echo "assembly 1615942301")

##############################################################################
# 15_YaHS_CuratedScaffolds.sh
##############################################################################
# after manual curation in juicebox, generate fixed scaffolded assembly

module load YaHS/1.1-GCC-11.3.0
module load Juicebox/2.20.00
module load SAMtools/1.16.1-GCC-11.3.0

juicer post -o fixedPAMFEO out_JBAT.review.FINISH.assembly out_JBAT.liftover.agp filteredcontigs.fa

##############################################################################
# 16_BUSCO_PostScaf.sh
##############################################################################
# busco

module load BUSCO/5.5.0-foss-2022a

busco -m genome -r -i PAMFEO_Assembly.fa -o buscoPostScaf -l insecta_odb10 --cpu 32 --download_path ~/Genome/busco_downloads