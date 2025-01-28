# Sections indicate separate cluster submission scripts

# Modules used:
BamTools/2.5.1
BLAST+/2.12.0
pigz/2.7
hifiasm/0.19.6
SeqKit/2.5.1
Kraken2/2.1.3
RepeatMasker/4.1.5
HISAT2/3n-20201216
SAMtools/1.16.1

# Folders
~/Genome/RawData # contains R1 and R2 Hi-C reads and HiFi reads
~/Genome/Transcriptome # contains R1 and R2 reads for 8 insects/3 tissues each from separate RNA-seq experiment

# Output
filteredcontigs.fa # primary assembly from hifiasm, filtered by PacBio coverage, RNA-seq coverage, and Kraken2 hits

##############################################################################
# 1_HiFiAdaptFilt.sh
##############################################################################
# removing adapters from ccs reads
# due to data transfer issues, I am starting partway through the process with the blocklist file generated by ARS-USDA

module load BamTools/2.5.1-GCC-10.2.0
module load BLAST+/2.12.0-gompi-2020b
module load pigz/2.7-GCCcore-10.2.0

cd ~/Genome/RawData

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

##############################################################################
# 2_HiFiASM.sh
##############################################################################
# running hifiasm with HiC data 

module load hifiasm/0.19.6-GCCcore-11.3.0

cd ~/Genome/hifiasm

hifiasm -o PAMFEO.asm -t 64 --h1 ~/Genome/RawData/HiCPAMFEO_R1.fastq.gz \
--h2 ~/Genome/RawData/HiCPAMFEO_R2.fastq.gz \
-f38 --primary ~/Genome/RawData/PAMFEO_allruns.fastq.gz

##############################################################################
# 3_GFA2Fasta_PBCovStat.sh
##############################################################################
# convert to fasta and get stats for HiFi read coverage

cd ~/Genome/hifiasm
mkdir hap

# Convert to fasta file
awk '/^S/{print ">"$2;print $3}' PAMFEO.asm.hic.p_ctg.gfa > PAMFEO_p.fa
awk '/^S/{print ">"$2;print $3}' PAMFEO.asm.hic.hap1.p_ctg.gfa > hap/PAMFEO_hap1.fa
awk '/^S/{print ">"$2;print $3}' PAMFEO.asm.hic.hap2.p_ctg.gfa > hap/PAMFEO_hap2.fa

# Obtain data for contig length and coverage
grep ^S PAMFEO.asm.hic.p_ctg.noseq.gfa > PBcov_primary.tsv
grep ^S PAMFEO.asm.hic.hap1.p_ctg.noseq.gfa > PBcov_hap1.tsv
grep ^S PAMFEO.asm.hic.hap2.p_ctg.noseq.gfa > PBcov_hap2.tsv

# evaluate contigs in Excel based on PacBio coverage, generate Figure 1 (edited in Illustrator for aesthetics)
# I kept contigs with coverage (rd:i:##) 6=< ## =<30

##############################################################################
# 4_LinearizeSplitFragment.sh
##############################################################################
# linearize fasta for easier grep-ing, split into contig files and fragment for Kraken2

module load SeqKit/2.5.1

cd ~/Genome/hifiasm
mkdir contig_ind # for individual contig .fa files
mkdir contig_frag # for fragmented contig .fa files

# Linearize genome so every contig sequence is one line long
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < PAMFEO_p.fa > int.fa
tail -n +2 int.fa > linear_PAMFEO.fa
rm int.fa

# filter by PacBio coverage
# Use excel function =TEXTJOIN("\|",TRUE,[range]) to combine contig names of ones you want to keep for use with grep
#	remember to use ' ' on set of contig names
#		grep pattern shortened for this script
grep -A 1 'contig1\|contig2\|contig3\|......' linear_PAMFEO.fa > filtcontigs.fa
grep -v '\-' filtcontigs.fa > PB_filtcontigs.fa # removes line between contigs (I could probably have used --no-group-separator flag)

# Place every contig into its own file
cd contig_ind
cat ~/Genome/hifiasm/PB_filtcontigs.fa| awk '{
        if (substr($0, 1, 1)==">") {filename=(substr($0,2) ".fa")}
        print $0 >> filename
        close(filename)
}'

# Make file of contig names for all .fa files in current directory
ls *.fa > files.txt 
cd -

# Use seqkit to split contigs into 100nt fragments that overlap by 50nt
for line in $(cat contig_ind/files.txt); do
seqkit sliding -s 50 -W 100 contig_ind/$line > contig_frag/frag"$line"
done

# Make file of fragmented contig names within contig_frag directory
cd contig_frag
ls *.fa > ~/Genome/hifiasm/fragfiles.txt

##############################################################################
# 5a_Kraken2.sh
##############################################################################
# run fragmented contig pseudo-reads through Kraken2
# run on cluster as array with array=1-(number_of_contigs):

module load Kraken2/2.1.3-gompi-2022a

cd ~/Genome/hifiasm

file=$(awk "NR==${SLURM_ARRAY_TASK_ID}" fragfiles.txt)
echo $file

kraken2 --use-names --threads 32 --db /db/kraken2 --report krakenreports/"$file"_report.txt \
--minimum-hit-groups 3 contig_frag/"$file" > krakenout/"$file".kraken

##############################################################################
# 5b_RepeatMasker_Simple.sh
##############################################################################
# mask intact contig files concurrently with Kraken2
# run as array with array=1-(number_of_contigs)

module load RepeatMasker/4.1.5-foss-2022a

cd ~/Genome/hifiasm
mkdir masked

file=$(awk "NR==${SLURM_ARRAY_TASK_ID}" contig_ind/files.txt)
echo $file

RepeatMasker -q -pa 32 -species insecta -dir masked contig_ind/$file

##############################################################################
# 6a_Hisat2BuildIndex.sh
##############################################################################
# prep for hisat

module load HISAT2/3n-20201216-gompi-2022a

cd ~/Genome/hifiasm

cat masked/*.fa.masked > maskedgenome.fa

hisat2-build maskedgenome.fa maskedgenome

##############################################################################
# 6b_Hisat2Array.sh
##############################################################################
# make array with RNA samples so each can be done separately
# RNA samples: 8x fat body (FB), 8x hindgut (HgT), 8x midgut (MgT) processed separately from this study

module load HISAT2/3n-20201216-gompi-2022a
module load SAMtools/1.16.1-GCC-11.3.0

cd ~/Genome/hifiasm

file=$(awk "NR==${SLURM_ARRAY_TASK_ID}" RNASamples.txt)

hisat2 -p 32 -x maskedgenome -k 3 --dta --no-unal \
-1 Transcriptome/PamRNAseqs/rep2"$file"_R1.fastq.gz -2 Transcriptome/PamRNAseqs/rep2"$file"_R2.fastq.gz \
--rg-id "$file" --rg SM:"$file" --summary-file hisat2summary_"$file".txt | \
samtools view -@32 -Sb - | \
samtools sort -@32 --write-index - -o filt_masked"$file"_sort.bam##idx##filt_masked"$file"_sort.bam.bai

##############################################################################
# 7a_mergeFB_alignments.sh
##############################################################################
# merge tissue-specific bam files and generate RNA coverage file

module load SAMtools/1.16.1-GCC-11.3.0

samtools merge -@ 16 -o merged_FB.bam \
filt_maskedXT1_FB1_sort.bam filt_maskedXT1_FB2_sort.bam filt_maskedXT1_FB3_sort.bam filt_maskedXT1_FB5_sort.bam \
filt_maskedXT5_FB1_sort.bam filt_maskedXT5_FB2_sort.bam filt_maskedXT5_FB3_sort.bam filt_maskedXT5_FB6_sort.bam

samtools coverage -o FBcov.tsv merged_FB.bam

##############################################################################
# 7b_mergeHgT_alignments.sh
##############################################################################
# merge tissue-specific bam files and generate RNA coverage file

module load SAMtools/1.16.1-GCC-11.3.0

samtools merge -@ 16 -o merged_HgT.bam \
filt_maskedXT1_HgT1_sort.bam filt_maskedXT1_HgT2_sort.bam filt_maskedXT1_HgT3_sort.bam filt_maskedXT1_HgT5_sort.bam \
filt_maskedXT5_HgT1_sort.bam filt_maskedXT5_HgT2_sort.bam filt_maskedXT5_HgT3_sort.bam filt_maskedXT5_HgT6_sort.bam

samtools coverage -o HgTcov.tsv merged_HgT.bam

##############################################################################
# 7c_mergeMgT_alignments.sh
##############################################################################
# merge tissue-specific bam files and generate RNA coverage file

module load SAMtools/1.16.1-GCC-11.3.0

samtools merge -@ 16 -o merged_MgT.bam \
filt_maskedXT1_MgT1_sort.bam filt_maskedXT1_MgT2_sort.bam filt_maskedXT1_MgT3_sort.bam filt_maskedXT1_MgT5_sort.bam \
filt_maskedXT5_MgT1_sort.bam filt_maskedXT5_MgT2_sort.bam filt_maskedXT5_MgT3_sort.bam filt_maskedXT5_MgT6_sort.bam

samtools coverage -o MgTcov.tsv merged_MgT.bam

##############################################################################
# 8_FilterandGrep.txt
##############################################################################
#In excel, mark which contigs you want to keep after masking repeats, running Kraken2, and running Hisat2 with RNA reads
# Use excel function =TEXTJOIN("\|",TRUE,[range]) to combine contig names of ones you want to keep for use with grep
#	remember to use ' ' on set of contig names
#		grep pattern shortened for this script

grep -A 1 'contig1\|contig2\|contig3\|......' linear_PAMFEO.fa > tempcontigs.fa
grep -v '\-' tempcontigs.fa > filteredcontigs.fa # removes line between contigs
rm tempcontigs.fa