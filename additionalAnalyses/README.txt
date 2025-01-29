These are the scripts used to perform additional analyses on the genome composition. 

TIDK.sh 
    telomere identification, ran on the initial scaffolded genome assembly and again after manual curating in Juicebox
RepeatModelMask.sh
    collection of bash scripts based on Daren Card's repeat modeling tutorial: https://darencard.net/blog/2022-07-09-genome-repeat-annotation/
    uses the scripts parseRM.pl, repclassifier.sh, and rmOutToGFF3custom.sh located in otherScripts folder
MitoHiFi.sh
    mitochondria identification from initial hifiasm assembly
