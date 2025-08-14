## open from terminal in orthofinder conda environment:
# conda activate orthofinder
# open -na rstudio

# libraries for GENESPACE
library("devtools")
library("BiocManager")
library("tidyverse")
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
#devtools::install_github("jtlovell/GENESPACE")
#BiocManager::install(c("Biostrings", "rtracklayer"))
library("GENESPACE")

# libraries for orthofinder and GO analysis
#BiocManager::install("clusterProfiler")
library("clusterProfiler")
library("ggprism")
library("data.table")

# Prep and run GENESPACE --------------------------------------------------
# NCBI format
parsedPaths <- parse_annotations(
  rawGenomeRepo = "~/genespace/genomeRepo",
  genomeDirs = c("BGer","CFor","CSec","DPun","ZooNev","Pamfeo"
                  ),
 #                 ,"PAm2022"),
  genomeIDs = c("BGer","CFor","CSec","DPun","ZooNev","Pamfeo"
                ), 
 #               ,"PAm2022"),
  gffString = "gff",
  faString = "fa",
  presets = "ncbi",
  genespaceWd = "~/genespace")

genomes2run <- c("BGer","CFor","CSec","DPun","ZooNev", "Pamfeo")
gpar <- init_genespace(
  wd = "~/genespace", 
  path2mcscanx = "~/genespace/MCScanX-master")
runGS <- run_genespace(gsParam = gpar) 

# GENESPACE Riparian Plots ------------------------------------------------
ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white"))
customPal <- colorRampPalette(
  c("darkorange", "skyblue", "darkblue", "purple", "darkred", "salmon"))
customPal2 <- colorRampPalette(
  c("#125A56", "#238F9D", "#C6DBED", "#F0E6B2","#FD9A44", "salmon"))

# comparing PAMFEO with available cockroach genomes/annotations
riparian_roaches <- plot_riparian(
  gsParam = runGS, 
  refGenome = "Pamfeo",
  genomeIDs = c("BGer","Pamfeo","DPun"),
  useRegions = TRUE,
  useOrder = FALSE,
  palette = customPal2,
  chrFill = "lightgrey",
  braidAlpha = .75,
  addThemes = ggthemes)

# comparing PAMFEO with available RefSeq termite genomes/annotations
riparian_termites <-plot_riparian(
  gsParam = runGS, 
  refGenome = "Pamfeo",
  genomeIDs = c("ZooNev","Pamfeo","CSec"),
  useRegions = TRUE,
  useOrder = FALSE,
  palette = customPal2,
  chrFill = "lightgrey",
  braidAlpha = .75,
  addThemes = ggthemes)

# comparison between this P. americana and the previous genome
# riparian_periplaneta <- plot_riparian(
#   gsParam = runGS, 
#   refGenome = "Pamfeo",
#   genomeIDs = c("PAm2022","Pamfeo"),
#   useRegions = TRUE,
#   useOrder = FALSE,
#   palette = customPal2,
#   chrFill = "lightgrey",
#   braidAlpha = .75,
#   addThemes = ggthemes)

# all Blattodea; not very informative without indicating particular genes of interest
ripDat2 <- plot_riparian(
  gsParam = runGS, 
  refGenome = "Pamfeo",
  genomeIDs = c("BGer","DPun","Pamfeo","CFor","CSec","ZooNev" ),
  useRegions = TRUE,
  useOrder = FALSE,
  palette = customPal2,
  chrFill = "lightgrey",
  braidAlpha = .75,
  addThemes = ggthemes)

# UpSetR ------------------------------------------------------------------
library(UpSetR)
genecount <- read.table("~/genespace/orthofinder/Results/Orthogroups/Orthogroups.GeneCount.tsv", header = TRUE, row.names = 1)
upset_gene <- genecount[,1:6] # removing "total" column
upset_gene[upset_gene > 0] <- 1 # change to presence/absence
upset(upset_gene, nsets = 8, sets = rev(c("ZooNev","CFor","CSec","Pamfeo", 
                                          # "PAm2022",  
                                          "BGer","DPun")), 
      keep.order = T, order.by = "freq", nintersects = 100)
write.csv(upset_gene, "upset_OGs.csv")

# attaching GO terms to OGs based on PAMFEO -------------------------------

# import modified GAF file for P. americana annotations
    # in Excel:
    # removed header
    # pulled out columns GeneID, Symbol, Qualifier, GO_ID, Aspect, and Gene_Name
    # deduplicated
pamfeoGAF <- read.csv("~/genespace/PAMFEO_NCBI_modifiedGAF.csv")

# import modified GFF file for P. americana gene descriptions
    # in Excel:
    # removed header
    # filtered for feature "gene"
    # split information associated with each gene by ";"
    # retained Symbol, chromosome placement, and description
    # deduplicated
pamfeoGFF <- read.csv("~/genespace/PAMFEO_NCBI_modifiedGFF.csv")

# import gene to orthogroup assignment
Orthogroups <- read.delim2("~/genespace/orthofinder/Results/Orthogroups/Orthogroups.tsv")

# identify OGs present in PAMFEO
allOG.PAMFEO <- rownames(upset_gene[which(upset_gene$Pamfeo == 1),]) # 10795 OGs

# pull out genes belonging to PAMFEO-containing OGs
allGenes.PAMFEO <- Orthogroups[which(Orthogroups$Orthogroup %in% allOG.PAMFEO),c("Orthogroup","Pamfeo")] %>% 
  remove_rownames() %>%
  column_to_rownames("Orthogroup") %>%
  apply(1,function(x) str_split(x,pattern = ", ")) %>%
  lapply(unlist) %>% 
  melt() 
colnames(allGenes.PAMFEO) <- c("Symbol", "Orthogroup")
length(unique(allGenes.PAMFEO$Symbol)) # 15926, no duplicates

# add gene description and other info to genes that are in OGs
allGeneswGFF.PAMFEO <- inner_join(allGenes.PAMFEO,pamfeoGFF)
setdiff(allGenes.PAMFEO$Symbol,allGeneswGFF.PAMFEO$Symbol) # noticed that the numbers did not match, caused by parentheses in GFF file
# after fixing parentheses, one mismatch between these files "beta_COP", but they both have a "betaCOP". Dropping. 
# 15925 genes

# combine GAF file with OG assignments and GO ID, adding GO term
allGAF.PAMFEO <- left_join(allGeneswGFF.PAMFEO,pamfeoGAF)
allTerm.PAMFEO <- go2term(allGAF.PAMFEO$GO_ID)
colnames(allTerm.PAMFEO) <- c("GO_ID","Term")
allGAF.wTerm <- full_join(allGAF.PAMFEO,allTerm.PAMFEO)

# save file for records
write.csv(allGAF.wTerm, "~/genespace/allGAFwTerm.csv")

# PAMFEO-unique OGs  ------------------------------------------------------

# identify OGs unique to PAMFEO and extract from modified GAF table
uPAMFEO.OGs <- rownames(upset_gene[which(rowSums(upset_gene)==1 & upset_gene$Pamfeo == 1),]) # 269 OGs
uPAMFEO.GAF <- allGAF.wTerm[which(allGAF.wTerm$Orthogroup %in% uPAMFEO.OGs),]
length(unique(uPAMFEO.GAF[!is.na(uPAMFEO.GAF$GO_ID),"Orthogroup"])) # 170 OGs have GO value

# save file for records
write.csv(uPAMFEO.GAF, "~/genespace/uniquePAMFEOgafOG.csv")

# weighting multi-GO term genes 
uPAMFEO.MultiTerm <- uPAMFEO.GAF[,c("Symbol","Aspect")] %>%
  dplyr::group_by_all() %>%
  count
uPAMFEO.MultiTerm$Weight <- 1/uPAMFEO.MultiTerm$n
uPAMFEO.weighted <- left_join(uPAMFEO.GAF,uPAMFEO.MultiTerm, by = c("Symbol","Aspect"))
uPAMFEO.Summarized <- uPAMFEO.weighted[,c("GO_ID","Aspect","Term","Weight")] %>%
  dplyr::group_by(pick(GO_ID,Aspect,Term)) %>%
  summarise(weightSum = sum(Weight))

# adding column to assist in plotting order
uPAMFEO.Summarized$AspectTerm <- paste(uPAMFEO.Summarized$Aspect,
                                       uPAMFEO.Summarized$Term, sep = ": ")
write.csv(uPAMFEO.Summarized, "~/genespace/summarizedGOTerms_unique.csv")

# remove NA row
uPAMFEO.Summarized <- uPAMFEO.Summarized[!is.na(uPAMFEO.Summarized$GO_ID),]

# exporting as PDF to adjust in Adobe Illustrator
pdf("~/genespace/draftUniqueOGGoTerms_Weight_14.pdf", width = 10, height = 10)
ggplot(uPAMFEO.Summarized[which(uPAMFEO.Summarized$weightSum > 14),], 
       mapping = aes(1,AspectTerm, size = weightSum, color = Aspect)) + 
  geom_point() + theme_prism()
dev.off()

# now using gene descriptions to capture genes without GO IDs
uPAMFEO.desc <- uPAMFEO.GAF[!duplicated(uPAMFEO.GAF[,c("Symbol","Orthogroup","description")]),c("Symbol","Orthogroup","description")]

# simplifying common descriptions, maybe better way to do this.
uPAMFEO.desc[grep("zinc finger",uPAMFEO.desc$description),"description"] <- "zinc finger"
uPAMFEO.desc[grep("uncharacterized",uPAMFEO.desc$description),"description"] <- "uncharacterized"
uPAMFEO.desc[grep("ankyrin",uPAMFEO.desc$description),"description"] <- "ankyrin"
uPAMFEO.desc[grep("^trypsin",uPAMFEO.desc$description),"description"] <- "trypsin"
uPAMFEO.desc[grep("odorant receptor",uPAMFEO.desc$description),"description"] <- "odorant receptor"
uPAMFEO.desc[grep("lipase",uPAMFEO.desc$description),"description"] <- "lipase"
uPAMFEO.desc[grep("very long chain fatty acid elongase",uPAMFEO.desc$description),"description"] <- "very long chain fatty acid elongase"
uPAMFEO.desc[grep("golgin subfamily",uPAMFEO.desc$description),"description"] <- "golgin subfamily"
uPAMFEO.desc[grep("chymotrypsin",uPAMFEO.desc$description),"description"] <- "chymotrypsin"
uPAMFEO.desc[grep("cuticle protein",uPAMFEO.desc$description),"description"] <- "cuticle protein"
uPAMFEO.desc[grep("cytochrome P450",uPAMFEO.desc$description),"description"] <- "cytochrome P450"

uPAMFEO.descSummary <- uPAMFEO.desc[,c("Orthogroup","description")] %>%
  dplyr::group_by_all() %>%
  count

uPAMFEO.descSummary2 <- uPAMFEO.descSummary %>%
  dplyr::group_by(description) %>%
  summarise(sum = sum(n))
plotThese <- uPAMFEO.descSummary2[which(uPAMFEO.descSummary2$sum >= 10),"description"]

# plotting to PDF, editing in Illustrator
pdf("~/genespace/uPAMFEOGeneNames.pdf")
ggplot(uPAMFEO.descSummary[which(uPAMFEO.descSummary$description %in% plotThese$description),],
       aes(Orthogroup,description,size = n)) + 
  geom_point()
dev.off()

# Shared OGs --------------------------------------------------------------
## using the GO annotations from the P. americana assembly to give broad overview of OG functions
# identify OGs shared by all groups
shared.OGs <- rownames(upset_gene[which(rowSums(upset_gene)==6 ),]) # 4532 OGs
shared.GAF <- allGAF.wTerm[which(allGAF.wTerm$Orthogroup %in% shared.OGs),]
length(unique(shared.GAF[!is.na(shared.GAF$GO_ID),"Orthogroup"])) # 3769 OGs have GO value

# save file for records
write.csv(shared.GAF, "~/genespace/sharedPAMFEOgafOG.csv")

# weighting multi-GO term genes 
shared.MultiTerm <- shared.GAF[,c("Symbol","Aspect")] %>%
  dplyr::group_by_all() %>%
  count
shared.MultiTerm$Weight <- 1/shared.MultiTerm$n
shared.weighted <- left_join(shared.GAF,shared.MultiTerm, by = c("Symbol","Aspect"))
shared.Summarized <- shared.weighted[,c("GO_ID","Aspect","Term","Weight")] %>%
  dplyr::group_by(pick(GO_ID,Aspect,Term)) %>%
  summarise(weightSum = sum(Weight))

# adding column to assist in plotting order
shared.Summarized$AspectTerm <- paste(shared.Summarized$Aspect,
                                      shared.Summarized$Term, sep = ": ")
write.csv(shared.Summarized, "~/genespace/summarizedGOTerms_shared.csv")

# remove NA row
shared.Summarized <- shared.Summarized[!is.na(shared.Summarized$GO_ID),]

# exporting as PDF to adjust in Adobe Illustrator
pdf("~/genespace/draftSharedOGGoTerms_Weight_14.pdf", width = 10, height = 10)
ggplot(shared.Summarized[which(shared.Summarized$weightSum > 14),], 
       mapping = aes(1,AspectTerm, size = weightSum, color = Aspect)) + 
  geom_point() + theme_prism()
dev.off()



# Session Info ------------------------------------------------------------
# R version 4.2.3 (2023-03-15)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS 14.0
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] GENESPACE_1.3.1       UpSetR_1.4.0          BiocGenerics_0.44.0   lubridate_1.9.3       forcats_1.0.0         stringr_1.5.1         dplyr_1.1.4           purrr_1.0.2           readr_2.1.5          
# [10] tidyr_1.3.1           tibble_3.2.1          ggplot2_3.5.1         tidyverse_2.0.0       BiocManager_1.30.25   devtools_2.4.5        usethis_3.1.0         data.table_1.15.4     ggprism_1.0.5        
# [19] clusterProfiler_4.6.2
# 
# loaded via a namespace (and not attached):
#   [1] circlize_0.4.16             shadowtext_0.1.4            fastmatch_1.1-4             plyr_1.8.9                  igraph_2.0.3                lazyeval_0.2.2              proxyC_0.4.1               
# [8] splines_4.2.3               BiocParallel_1.32.6         GenomeInfoDb_1.34.9         digest_0.6.35               foreach_1.5.2               yulab.utils_0.2.0           htmltools_0.5.8.1          
# [15] GOSemSim_2.24.0             viridis_0.6.5               GO.db_3.16.0                magrittr_2.0.3              memoise_2.0.1               tm_0.7-13                   cluster_2.1.6              
# [22] doParallel_1.0.17           tzdb_0.4.0                  remotes_2.5.0               ComplexHeatmap_2.14.0       Biostrings_2.66.0           annotate_1.76.0             graphlayouts_1.1.1         
# [29] matrixStats_1.3.0           R.utils_2.12.3              timechange_0.3.0            enrichplot_1.18.4           colorspace_2.1-0            blob_1.2.4                  ggrepel_0.9.5              
# [36] crayon_1.5.3                RCurl_1.98-1.14             jsonlite_1.8.8              ggsankeyfier_0.1.8          scatterpie_0.2.4            iterators_1.0.14            ape_5.8                    
# [43] glue_1.7.0                  polyclip_1.10-6             gtable_0.3.6                zlibbioc_1.44.0             XVector_0.38.0              GetoptLong_1.0.5            DelayedArray_0.24.0        
# [50] pkgbuild_1.4.6              shape_1.4.6.1               scales_1.3.0                DOSE_3.24.2                 DBI_1.2.3                   miniUI_0.1.1.1              Rcpp_1.0.12                
# [57] viridisLite_0.4.2           xtable_1.8-4                clue_0.3-65                 gridGraphics_0.5-1          tidytree_0.4.6              bit_4.0.5                   stats4_4.2.3               
# [64] profvis_0.3.8               htmlwidgets_1.6.4           httr_1.4.7                  fgsea_1.24.0                RColorBrewer_1.1-3          ellipsis_0.3.2              R.methodsS3_1.8.2          
# [71] urlchecker_1.0.1            pkgconfig_2.0.3             XML_3.99-0.16.1             farver_2.1.2                utf8_1.2.4                  locfit_1.5-9.9              labeling_0.4.3             
# [78] ggplotify_0.1.2             tidyselect_1.2.1            rlang_1.1.3                 reshape2_1.4.4              later_1.3.2                 AnnotationDbi_1.60.2        munsell_0.5.1              
# [85] tools_4.2.3                 cachem_1.1.0                downloader_0.4              cli_3.6.2                   dbscan_1.1-12               generics_0.1.3              RSQLite_2.3.7              
# [92] gson_0.1.0                  fastmap_1.2.0               ggtree_3.6.2                org.Hs.eg.db_3.16.0         bit64_4.0.5                 fs_1.6.4                    tidygraph_1.3.1            
# [99] KEGGREST_1.38.0             ggraph_2.2.1                nlme_3.1-164                mime_0.12                   slam_0.1-50                 R.oo_1.27.0                 aplot_0.2.4                
# [106] xml2_1.3.6                  compiler_4.2.3              rstudioapi_0.17.1           png_0.1-8                   treeio_1.22.0               tweenr_2.0.3                geneplotter_1.76.0         
# [113] stringi_1.8.4               lattice_0.22-6              Matrix_1.6-5                vctrs_0.6.5                 pillar_1.10.1               lifecycle_1.0.4             GlobalOptions_0.1.2        
# [120] cowplot_1.1.3               bitops_1.0-7                httpuv_1.6.15               patchwork_1.3.0             GenomicRanges_1.50.2        qvalue_2.30.0               R6_2.5.1                   
# [127] promises_1.3.0              gridExtra_2.3               IRanges_2.32.0              sessioninfo_1.2.3           codetools_0.2-20            MASS_7.3-60.0.1             pkgload_1.4.0              
# [134] SummarizedExperiment_1.28.0 rjson_0.2.21                DESeq2_1.38.3               withr_3.0.2                 S4Vectors_0.36.2            GenomeInfoDbData_1.2.9      parallel_4.2.3             
# [141] hms_1.1.3                   ggfun_0.1.8                 HDO.db_0.99.1               MatrixGenerics_1.10.0       ggforce_0.4.2               NLP_0.3-2                   Biobase_2.58.0             
# [148] shiny_1.9.1                
