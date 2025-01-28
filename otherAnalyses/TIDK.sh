# get program here: https://github.com/tolkit/telomeric-identifier

conda install -c bioconda tidk
conda activate tidk

tidk explore --minimum 5 --maximum 12 PAMFEO_Assembly.fa > tidk_explore.tsv

# chose 5 promising options to find and plot
tidk search --extension tsv --string AACCTAACCT --output AACCTAACCT --dir telo PAMFEO_Assembly.fa
tidk search --extension tsv --string TTAGG --output TTAGG --dir telo PAMFEO_Assembly.fa
tidk search --extension tsv --string AAAATTC --output AAAATTC --dir telo PAMFEO_Assembly.fa
tidk search --extension tsv --string AACTCAGCG --output AACTCAGCG --dir telo PAMFEO_Assembly.fa
tidk search --extension tsv --string AATACAATAC --output AATACAATAC --dir telo PAMFEO_Assembly.fa

tidk plot --tsv telo/AACCTAACCT_telomeric_repeat_windows.tsv --output AACCTAACCT_Final
tidk plot --tsv telo/TTAGG_telomeric_repeat_windows.tsv --output TTAGG_Final
tidk plot --tsv telo/AAAATTC_telomeric_repeat_windows.tsv --output AAAATTC_Final
tidk plot --tsv telo/AACTCAGCG_telomeric_repeat_windows.tsv --output AACTCAGCG_Final
tidk plot --tsv telo/AATACAATAC_telomeric_repeat_windows.tsv --output AATACAATAC_Final