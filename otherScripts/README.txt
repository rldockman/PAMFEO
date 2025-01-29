These are additional scripts I used directly and without modification from other sources. 
Arima Genomics HiC mapping pipeline: https://github.com/ArimaGenomics/mapping_pipeline
  filter_five_end.pl
      after aligning HiC reads to the genome, this filters out chimeric alignments aligned in the 3' orientation (since that indicates a ligation junction, not a contiguous read)
  two_read_bam_combiner.pl 
      combine filtered forward and reverse mapped HiC reads to a single bam file
  get_stats.pl
      get stats on inter- and intra-contig read pairs remaining after filtering and duplicate removal
