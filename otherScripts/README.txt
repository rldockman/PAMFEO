These are additional scripts I used directly and without modification from other sources. 

Arima Genomics HiC mapping pipeline: https://github.com/ArimaGenomics/mapping_pipeline
  filter_five_end.pl
      after aligning HiC reads to the genome, this filters out chimeric alignments aligned in the 3' orientation (since that indicates a ligation junction, not a contiguous read)
  two_read_bam_combiner.pl 
      combine filtered forward and reverse mapped HiC reads to a single bam file
  get_stats.pl
      get stats on inter- and intra-contig read pairs remaining after filtering and duplicate removal

Daren Card's guide to repeat modeling and masking: https://darencard.net/
    repclassifer.sh
        ran iteratively to classify known and unknown repeats within the genome
    rmOutToGFF3custom.sh
        converts .out file to GFF3, allows for better repeat family summarization of RepeatModeler REs

Aurelie Kapusta's RepeatMasker parsing script: https://github.com/4ureliek
    parseRM.pl
        parses through repeat results to generate repeat landscape files
