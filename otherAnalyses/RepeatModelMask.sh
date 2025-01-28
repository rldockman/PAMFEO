# based off of: https://darencard.net/blog/2022-07-09-genome-repeat-annotation/

##############################################################################
# Model and Classify Repeats
##############################################################################
# used interactive node for most steps

ml RepeatModeler/2.0.4-foss-2022a
ml RepeatMasker/4.1.4-foss-2022a
ml bioawk/1.0-GCC-11.2.0
ml SeqKit/2.5.1

# repeat modeler
BuildDatabase -name Periplaneta_americana -engine ncbi PAMFEO1_hifiasm_priV1.fa
RepeatModeler -pa 32 -engine ncbi -database Periplaneta_americana 2>&1 | tee 00_repeatmodeler.log

# separate modeled repeat families into "known" and "unknown"
cat Periplaneta_americana-families.fa | seqkit fx2tab | awk '{ print "PAMFEO_"$0 }' | seqkit tab2fx > Periplaneta_americana-families.PAMFEO.fa
cat Periplaneta_americana-families.PAMFEO.fa | seqkit fx2tab | grep -v "Unknown" | seqkit tab2fx > Periplaneta_americana-families.PAMFEO.fa.known
cat Periplaneta_americana-families.PAMFEO.fa | seqkit fx2tab | grep "Unknown" | seqkit tab2fx > Periplaneta_americana-families.PAMFEO.fa.unknown

# quantify number of classified elements
grep -c ">" Periplaneta_americana-families.PAMFEO.fa.known # 1492

# quantify number of unknown elements
grep -c ">" Periplaneta_americana-families.PAMFEO.fa.unknown # 2954

# classifying unknowns (-u): run with 3 threads/cores (-t) and using the Insecta elements (-d) from Repbase/Dfam 
# 	and known elements (-k) from the same reference genome; append newly identified elements to the existing known 
# 	element library (-a) and write results to an output directory (-o)
bash repclassifier.sh -t 8 -d Insecta -u Periplaneta_americana-families.PAMFEO.fa.unknown \
-k Periplaneta_americana-families.PAMFEO.fa.known -a Periplaneta_americana-families.PAMFEO.fa.known \
-o round-1_RepbaseInsecta-Self

# classifying unknowns (-u): run with 3 threads/cores (-t) and using only the known elements (-k) from the 
# 	same reference genome; append newly identified elements to the existing known element library (-a) and 
#	write results to an output directory (-o). No Repbase classification is used here.
bash repclassifier.sh -t 8 -u round-1_RepbaseInsecta-Self/round-1_RepbaseInsecta-Self.unknown \
-k round-1_RepbaseInsecta-Self/round-1_RepbaseInsecta-Self.known \
-a round-1_RepbaseInsecta-Self/round-1_RepbaseInsecta-Self.known -o round-2_Self

# repeat until no more unknown elements are classified

bash repclassifier.sh -t 8 -u RepeatModelMask/round-2_Self/round-2_Self.unknown \
-k RepeatModelMask/round-2_Self/round-2_Self.known \
-a RepeatModelMask/round-2_Self/round-2_Self.known -o round-3_Self

bash repclassifier.sh -t 8 -u round-3_Self/round-3_Self.unknown \
-k round-3_Self/round-3_Self.known \
-a round-3_Self/round-3_Self.known -o round-4_Self

bash repclassifier.sh -t 8 -u round-4_Self/round-4_Self.unknown \
-k round-4_Self/round-4_Self.known \
-a round-4_Self/round-4_Self.known -o round-5_Self

bash repclassifier.sh -t 8 -u round-5_Self/round-5_Self.unknown \
-k round-5_Self/round-5_Self.known \
-a round-5_Self/round-5_Self.known -o round-6_Self

bash repclassifier.sh -t 8 -u round-6_Self/round-6_Self.unknown \
-k round-6_Self/round-6_Self.known \
-a round-6_Self/round-6_Self.known -o round-7_Self

# round_7 produced no new families, so stopped here

mkdir -p logs 01_simple_out 02_Insecta_out 03_known_out 04_unknown_out

cp PAMFEO1_hifiasm_priV1.fa Periplaneta_americana.fasta

##############################################################################
# MultiRound_RM.sh
##############################################################################
# run through rounds of repeatmasker

ml SeqKit/2.5.1
ml RepeatModeler/2.0.4-foss-2022a
ml RepeatMasker/4.1.4-foss-2022a
ml bioawk/1.0-GCC-11.2.0

# round 1: annotate/mask simple repeats
RepeatMasker -pa 32 -a -e ncbi -dir 01_simple_out -noint -xsmall Periplaneta_americana.fasta 2>&1 | tee logs/01_simplemask.log
#rename
rename fasta simple_mask 01_simple_out/Periplaneta_americana*
rename .masked .masked.fasta 01_simple_out/Periplaneta_americana*

# round 2: annotate/mask Insecta elements sourced from Repbase using output from 1st round of RepeatMasker
RepeatMasker -pa 32 -a -e ncbi -dir 02_Insecta_out -nolow \
-species Insecta 01_simple_out/Periplaneta_americana.simple_mask.masked.fasta 2>&1 | tee logs/02_Insectamask.log

# round 2: rename outputs
rename simple_mask.masked.fasta Insecta_mask 02_Insecta_out/Periplaneta_americana*
rename .masked .masked.fasta 02_Insecta_out/Periplaneta_americana*

# round 3: annotate/mask known elements sourced from species-specific de novo repeat library using output from 2nd round of RepeatMasker
RepeatMasker -pa 32 -a -e ncbi -dir 03_known_out -nolow \
-lib round-7_Self/round-7_Self.known \
02_Insecta_out/Periplaneta_americana.Insecta_mask.masked.fasta 2>&1 | tee logs/03_knownmask.log

# round 3: rename outputs
rename Insecta_mask.masked.fasta known_mask 03_known_out/Periplaneta_americana*
rename .masked .masked.fasta 03_known_out/Periplaneta_americana*

# round 4: annotate/mask unknown elements sourced from species-specific de novo repeat library using output from 3nd round of RepeatMasker
RepeatMasker -pa 32 -a -e ncbi -dir 04_unknown_out -nolow \
-lib round-7_Self/round-7_Self.unknown \
03_known_out/Periplaneta_americana.known_mask.masked.fasta 2>&1 | tee logs/04_unknownmask.log

# round 4: rename outputs
rename known_mask.masked.fasta unknown_mask 04_unknown_out/Periplaneta_americana*
rename .masked .masked.fasta 04_unknown_out/Periplaneta_americana*

# create directory for full results
mkdir -p 05_full_out

# combine full RepeatMasker result files - .cat.gz
cat 01_simple_out/Periplaneta_americana.simple_mask.cat.gz \
02_Insecta_out/Periplaneta_americana.Insecta_mask.cat.gz \
03_known_out/Periplaneta_americana.known_mask.cat.gz \
04_unknown_out/Periplaneta_americana.unknown_mask.cat.gz \
> 05_full_out/Periplaneta_americana.full_mask.cat.gz

# combine RepeatMasker tabular files for all repeats - .out
cat 01_simple_out/Periplaneta_americana.simple_mask.out \
<(cat 02_Insecta_out/Periplaneta_americana.Insecta_mask.out | tail -n +4) \
<(cat 03_known_out/Periplaneta_americana.known_mask.out | tail -n +4) \
<(cat 04_unknown_out/Periplaneta_americana.unknown_mask.out | tail -n +4) \
> 05_full_out/Periplaneta_americana.full_mask.out

# copy RepeatMasker tabular files for simple repeats - .out
cat 01_simple_out/Periplaneta_americana.simple_mask.out > 05_full_out/Periplaneta_americana.simple_mask.out

# combine RepeatMasker tabular files for complex, interspersed repeats - .out
cat 02_Insecta_out/Periplaneta_americana.Insecta_mask.out \
<(cat 03_known_out/Periplaneta_americana.known_mask.out | tail -n +4) \
<(cat 04_unknown_out/Periplaneta_americana.unknown_mask.out | tail -n +4) \
> 05_full_out/Periplaneta_americana.complex_mask.out

# combine RepeatMasker repeat alignments for all repeats - .align
cat 01_simple_out/Periplaneta_americana.simple_mask.align \
02_Insecta_out/Periplaneta_americana.Insecta_mask.align \
03_known_out/Periplaneta_americana.known_mask.align \
04_unknown_out/Periplaneta_americana.unknown_mask.align \
> 05_full_out/Periplaneta_americana.full_mask.align

# resummarize repeat compositions from combined analysis of all RepeatMasker rounds
ProcessRepeats -a -species Insecta 05_full_out/Periplaneta_americana.full_mask.cat.gz 2>&1 | tee logs/05_fullmask.log

##############################################################################
# SoftMask.sh
##############################################################################
# soft masking the full repeat file for annotation

ml SeqKit/2.5.1
ml RepeatMasker/4.1.5-foss-2022a
ml BEDTools/2.30.0-GCC-11.3.0

bash rmOutToGFF3custom.sh -o RepeatModelMask/05_full_out/Periplaneta_americana.full_mask.out > FullRepeatMasking.gff3

bedtools maskfasta -soft -fi Periplaneta_americana.fasta -bed 05_full_out/Periplaneta_americana.simple_mask.gff3 \
-fo 05_full_out/Periplaneta_americana.fullmask.soft.fasta