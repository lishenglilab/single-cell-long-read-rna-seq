library(future)
library(IsoformSwitchAnalyzeR)
library(tidyverse)
library(muscat)
library(Seurat)
library(DESeq2)
library(BSgenome.Hsapiens.UCSC.hg38)
options(future.globals.maxSize = 80000000 * 1024^2)
plan(multisession, workers = 15)
sce.nano = read_rds('Result/AD_Keratinocytes_nanopore.rds')

if (T) {
  gtf <- rtracklayer::import('Reference/isoform_annotated.merged.sorted.gtf')
  gtf_df <- as.data.frame(gtf)
  gtf_df$transcript_id = gsub('_','-',gtf_df$transcript_id)

  filtered_gtf_df <- gtf_df[gtf_df$transcript_id %in% rownames(sce.nano), ]
  

  filtered_gtf <- makeGRangesFromDataFrame(filtered_gtf_df, keep.extra.columns = TRUE)

  rtracklayer::export(filtered_gtf, 'Reference/filtered_isoform_annotated.merged.sorted.gtf')
}


DefaultAssay(sce.nano) = 'RNA'
sce.nano <- as.SingleCellExperiment(sce.nano)
sce.nano <- prepSCE(sce.nano,
                     kid = "k_sub_ill",
                     gid = "stage",
                     sid = "orig.ident",
                     drop=T)

pb <- aggregateData(sce.nano,
                    assay = "counts", fun = "sum",
                    by = c("sample_id"))




isoformCountMatrix = as.data.frame(pb@assays@data@listData[[1]])
isoformCountMatrix = rownames_to_column(isoformCountMatrix,var = 'isoform_id')
isoformExonAnnoation =  'Reference/filtered_isoform_annotated.merged.sorted.gtf'
isoformNtFasta = 'Reference/transcript_assembly.fa'
localSwitchList <- importGTF(pathToGTF = isoformExonAnnoation, 
                                              addAnnotatedORFs = FALSE, removeTECgenes = FALSE, 
                                              quiet = TRUE)
isoAnnot <- unique(localSwitchList$isoformFeatures[, 
                                                   c("gene_id", "isoform_id", "gene_name")])
isoformCountMatrix = isoformCountMatrix[isoformCountMatrix$isoform_id %in% isoAnnot$isoform_id,]
geneCountMatrix = isoformToGeneExp(isoformCountMatrix,localSwitchList)
geneCountMatrix = select(geneCountMatrix,'gene_id',everything())

# DE analysis
run_deseq2 <- function(expression_matrix, sample_groups,ref) {
  rownames(expression_matrix) = expression_matrix[,1]
  expression_matrix[,1] = NULL
  if (length(sample_groups) != ncol(expression_matrix)) {
    stop("样本分组信息的长度必须与表达矩阵的列数一致")
  }
  
  sample_info <- data.frame(
    row.names = colnames(expression_matrix),
    condition = relevel(factor(sample_groups),ref = ref)
  )
  
  dds <- DESeqDataSetFromMatrix(countData = expression_matrix,
                                colData = sample_info,
                                design = ~ condition)

  dds <- DESeq(dds)

  results <- as.data.frame(results(dds))
  results$padj[is.na(results$padj)] = 1
  return(results)
}
sample_groups = c(rep(c('AD','Normal'),each = 3))
DEG = run_deseq2(geneCountMatrix,sample_groups,'Normal')
DET = run_deseq2(isoformCountMatrix,sample_groups,'Normal')

designMatrix <- data.frame(
  sampleID = colnames(isoformCountMatrix)[-1],
  condition = c(rep(c('AD','Normal'),each = 3))
)
comparisonsToMake <- data.frame(
  condition_1 = 'Normal',
  condition_2 = 'AD'
)

aSwitchList <- importRdata(
  isoformCountMatrix = isoformCountMatrix,
  designMatrix = designMatrix,
  isoformExonAnnoation = isoformExonAnnoation,
  comparisonsToMake = comparisonsToMake
)
### Filter
aSwitchList <- preFilter( aSwitchList )
### Test for isoform switches
aSwitchList <- isoformSwitchTestDEXSeq( aSwitchList )
extractSwitchSummary(aSwitchList)
#Part2
aSwitchList = analyzeORF(aSwitchList,BSgenome.Hsapiens.UCSC.hg38)
aSwitchList  = extractSequence(aSwitchList,genomeObject = BSgenome.Hsapiens.UCSC.hg38,pathToOutput  = 'Result/Consequence/',
                               outputPrefix = 'ALL_AD_Normal',onlySwitchingGenes=FALSE)

### Add annotation
aSwitchList <- analyzeCPC2( aSwitchList,
                            pathToCPC2resultFile = "Result/Consequence/extdata/result_cpc2.txt",removeNoncodinORFs = F)

aSwitchList <- analyzePFAM( aSwitchList,
                            pathToPFAMresultFile = "Result/Consequence/extdata/pfam_results_97.txt")
aSwitchList <- analyzeSignalP( aSwitchList,
                               pathToSignalPresultFile = "Result/Consequence/extdata/signalP_results_summary.signalp5")

aSwitchList <- analyzeIUPred2A(aSwitchList,
                               pathToIUPred2AresultFile = "Result/Consequence/extdata/iupred2a_result.txt")


### Analyse consequences
aSwitchList = analyzeIntronRetention(aSwitchList,onlySwitchingGenes = F)
saveRDS(aSwitchList,'Result/aSwitchList.rds')

df = read.table('Result/dominant_transcripts_pct_ad_normal.txt', sep = '\t', header = TRUE)
df$nor_pct[is.na(df$nor_pct)] = 0
df$ad_pct[is.na(df$ad_pct)] = 0
df$lgFC = log2((df$ad_pct + 0.000001) / (df$nor_pct + 0.000001))
df = na.omit(df)
res = df[abs(df$lgFC) >= 1,]
res$lgFC = log2(res$ad_pct / res$nor_pct)
res$change = ifelse(res$lgFC > 0, 'Up', 'Down')

# Generate pairwise comparisons
res_conse = res %>% 
  group_by(gene_id) %>% 
  filter(any(change == 'Up') & any(change == 'Down'))

pairwiseIsoComparison = res_conse %>%
  dplyr::group_by(gene_id) %>%
  dplyr::summarize(
    isoformUpregulated = list(transcript_id[change == "Up"]),
    isoformDownregulated = list(transcript_id[change == "Down"])
  ) %>%
  unnest(isoformUpregulated) %>%
  unnest(isoformDownregulated) %>%
  ungroup()

pairwiseIsoComparison$comparison = 1:nrow(pairwiseIsoComparison)
consequencesToAnalyze = c('tss','tts','last_exon','isoform_length','exon_number','intron_structure','ORF_length', '5_utr_seq_similarity', '5_utr_length', '3_utr_seq_similarity', '3_utr_length','coding_potential','ORF_seq_similarity','NMD_status','domains_identified','signal_peptide_identified')
# 
pairwiseIsoComparisonUniq <- unique(pairwiseIsoComparison[, c("isoformUpregulated", "isoformDownregulated")])
#

pairwiseIsoComparisonUniq$comparison <- 1:nrow(pairwiseIsoComparisonUniq)
minimumSwitchList <- makeMinimumSwitchList(orgSwitchList = aSwitchList, 
                                           isoformsToKeep = unique(c(pairwiseIsoComparisonUniq$isoformUpregulated, 
                                                                     pairwiseIsoComparisonUniq$isoformDownregulated)))

pairwiseIsoComparisonUniq = pairwiseIsoComparisonUniq[pairwiseIsoComparisonUniq$isoformUpregulated %in% minimumSwitchList$isoformFeatures$isoform_id & pairwiseIsoComparisonUniq$isoformDownregulated %in% minimumSwitchList$isoformFeatures$isoform_id,]

# Set cutoffs
ntCutoff = 50
ntFracCutoff = NULL
ntJCsimCutoff = 0.8
AaCutoff = 10
AaFracCutoff = 0.8
AaJCsimCutoff = 0.9

# Initialize parallel backend
numCores <- 27
cl <- makeCluster(numCores)
registerDoParallel(cl)

# Set up progress handling
handlers("progress")

# Run the function in parallel with progress tracking
with_progress({
  p <- progressor(along = pairwiseIsoComparisonUniq$comparison)
  consequencesOfIsoformSwitching <- foreach(aDF = iter(pairwiseIsoComparisonUniq, by = "row"), .packages = c("IsoformSwitchAnalyzeR"), .combine = 'rbind', .inorder = FALSE, .multicombine = TRUE, .export = c("compareAnnotationOfTwoIsoforms", "minimumSwitchList", "consequencesToAnalyze", "ntCutoff", "ntFracCutoff", "ntJCsimCutoff", "AaCutoff", "AaFracCutoff", "AaJCsimCutoff")) %dopar% {
    p() # Update progress
    compareAnnotationOfTwoIsoforms(
      switchAnalyzeRlist = minimumSwitchList,
      consequencesToAnalyze = consequencesToAnalyze,
      upIso = aDF$isoformUpregulated,
      downIso = aDF$isoformDownregulated,
      ntCutoff = ntCutoff,
      ntFracCutoff = ntFracCutoff,
      ntJCsimCutoff = ntJCsimCutoff,
      AaCutoff = AaCutoff,
      AaFracCutoff = AaFracCutoff,
      AaJCsimCutoff = AaJCsimCutoff,
      addDescription = TRUE,
      testInput = FALSE
    )
  }
})

# Stop the cluster
stopCluster(cl)
consequencesOfIsoformSwitching = left_join(consequencesOfIsoformSwitching,id_map[1:3],by = c('isoformDownregulated'='transcript_id'))
newOrder <- na.omit(match(c("gene_ref", "gene_id", "gene_name", 
                            "condition_1", "condition_2", "isoformUpregulated", 
                            "isoformDownregulated", "iso_ref_up", "iso_ref_down", 
                            "featureCompared", "isoformsDifferent", "switchConsequence"), 
                          colnames(consequencesOfIsoformSwitching)))
consequencesOfIsoformSwitchingDfcomplete <- consequencesOfIsoformSwitching[, newOrder]
consequencesOfIsoformSwitchingDfcomplete <- consequencesOfIsoformSwitchingDfcomplete[order(consequencesOfIsoformSwitchingDfcomplete$gene_id,consequencesOfIsoformSwitchingDfcomplete$isoformUpregulated, consequencesOfIsoformSwitchingDfcomplete$isoformDownregulated), 
]
write.table(consequencesOfIsoformSwitchingDfcomplete,'Result/consequencesOfIsoformSwitching.txt',sep = '\t',row.names = F,quote = F)
plot_data = as.data.frame(table(consequencesOfIsoformSwitchingDfcomplete$switchConsequence))
colnames(plot_data) = c('Functional_consequences','Number')
ggplot(plot_data,aes(Functional_consequences,Number))+
  geom_bar(stat = 'identity')+
  geom_text(aes(Functional_consequences,Number,label = Number),vjust = -0.3,size = 2.5)+
  theme_bw()+
  scale_y_continuous(expand = c(0,0),limits = c(0, max(plot_data$Number*1.1)))+
  theme(axis.text.x = element_text(angle = 45,vjust = 0.6),
        plot.margin = margin(t = 10, r = 15, b = 10, l = 10))
ggsave('Result/Functional_consequences_bar.pdf',width = 10,height = 6)

consequencesAnalyzed <- unique(consequencesOfIsoformSwitchingDfcomplete$featureCompared)
consequencesToAnalyze <- consequencesAnalyzed
localConseq <- consequencesOfIsoformSwitchingDfcomplete[which(!is.na(consequencesOfIsoformSwitchingDfcomplete$switchConsequence)), 
]
localConseq <- localConseq[which(!grepl("switch", localConseq$switchConsequence)), 
]
levelList <- list(tss = c("Tss more upstream", "Tss more downstream"), 
                  tts = c("Tts more downstream", "Tts more upstream"), 
                  last_exon = c("Last exon more downstream", "Last exon more upstream"), 
                  isoform_length = c("Length gain", "Length loss"), 
                  isoform_seq_similarity = c("Length gain", "Length loss"), 
                  exon_number = c("Exon gain", "Exon loss"), intron_retention = c("Intron retention gain", 
                                                                                  "Intron retention loss"), ORF_length = c("ORF is longer", 
                                                                                                                           "ORF is shorter"), ORF = c("Complete ORF loss", 
                                                                                                                                                      "Complete ORF gain"), x5_utr_length = c("5UTR is longer", 
                                                                                                                                                                                              "5UTR is shorter"), x3_utr_length = c("3UTR is longer", 
                                                                                                                                                                                                                                    "3UTR is shorter"), NMD_status = c("NMD sensitive", 
                                                                                                                                                                                                                                                                       "NMD insensitive"), coding_potential = c("Transcript is coding", 
                                                                                                                                                                                                                                                                                                                "Transcript is Noncoding"), domains_identified = c("Domain gain", 
                                                                                                                                                                                                                                                                                                                                                                   "Domain loss"), domain_length = c("Domain length gain", 
                                                                                                                                                                                                                                                                                                                                                                                                     "Domain length loss"), domain_isotype = c("Domain non-reference isotype gain", 
                                                                                                                                                                                                                                                                                                                                                                                                                                               "Domain non-reference isotype loss"), IDR_identified = c("IDR gain", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        "IDR loss"), IDR_length = c("IDR length gain", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    "IDR length loss"), IDR_type = c("IDR w binding region gain", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     "IDR w binding region loss"), signal_peptide_identified = c("Signal peptide gain", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 "Signal peptide loss"), sub_cell_location = c("SubCell location gain", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               "SubCell location loss"), sub_cell_shift_to_cell_membrane = c("SubCell location memb gain", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             "SubCell location memb loss"), sub_cell_shift_to_cytoplasm = c("SubCell location cyto gain", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            "SubCell location cyto loss"), sub_cell_shift_to_nucleus = c("SubCell location nucl gain", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         "SubCell location nucl loss"), sub_cell_shift_to_Extracellular = c("SubCell location ext cell gain", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            "SubCell location ext cell loss"), isoform_topology = c("Topology complexity gain", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    "Topology complexity loss"), extracellular_region_count = c("Extracellular region gain", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                "Extracellular region loss"), intracellular_region_count = c("Intracellular region gain", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             "Intracellular region loss"), extracellular_region_length = c("Extracellular length gain", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           "Extracellular length loss"), intracellular_region_length = c("Intracellular length gain", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         "Intracellular length loss"))
levelListDf <- plyr::ldply(levelList, function(x) data.frame(feature = x, 
                                                             stringsAsFactors = FALSE))
localConseq$conseqPair <- levelListDf$.id[match(localConseq$switchConsequence, 
                                                levelListDf$feature)]
localConseq <- localConseq[which(!is.na(localConseq$conseqPair)), 
]
localConseq <- localConseq[which(localConseq$featureCompared %in% 
                                   consequencesToAnalyze), ]
aSwitchList = read_rds('Result/aSwitchList.rds')
sigIso = aSwitchList$isoformFeatures[aSwitchList$isoformFeatures$isoform_id %in%unique(c( consequencesOfIsoformSwitchingDfcomplete$isoformUpregulated,consequencesOfIsoformSwitchingDfcomplete$isoformDownregulated)),c("iso_ref", "gene_ref")]
localConseq$condition_1 = 'Normal'
localConseq$condition_2 = 'AD'
consequenceBalance <- plyr::ddply(.data = localConseq, .variables = c("condition_1", "condition_2", "conseqPair"), .fun = function(aDF) {
  localLvl <- sort(levelList[[aDF$conseqPair[1]]])
  aDF$switchConsequence <- factor(aDF$switchConsequence,levels = localLvl)
  df2 <- aDF[which(!is.na(aDF$switchConsequence)),
             c("gene_id", "switchConsequence")]
  localNumber <- plyr::ddply(df2, .drop = FALSE, .variables = "switchConsequence",function(x) {
    data.frame(Freq = length(unique(x$gene_id)))
    })
  colnames(localNumber)[1] <- "Var1"
  if (nrow(localNumber) == 2) {
    localTest <- suppressWarnings(stats::binom.test(localNumber$Freq[1],
                                                    sum(localNumber$Freq)))
    localRes <- data.frame(feature = stringr::str_c(localNumber$Var1[1],
                                                    " (paired with ",
                                                    localNumber$Var1[2], ")"),
                           propOfRelevantEvents = localTest$estimate, stringsAsFactors = FALSE)
    localRes$propCiLo <- min(localTest$conf.int)
    localRes$propCiHi <- max(localTest$conf.int)
    localRes$propPval <- localTest$p.value
    }else {
      warning("Somthing strange happend - contact developer with reproducible example")
      }
  localRes$nUp <- localNumber$Freq[which(localNumber$Var1 == levels(localNumber$Var1)[1])]
  localRes$nDown <- localNumber$Freq[which(localNumber$Var1 == levels(localNumber$Var1)[2])]
  return(localRes)
  })
consequenceBalance$propQval <- p.adjust(consequenceBalance$propPval, 
                                        method = "fdr")
consequenceBalance$Significant <- consequenceBalance$propQval < 0.05
consequenceBalance$Significant <- factor(consequenceBalance$Significant, 
                                         levels = c(FALSE, TRUE))
xText <- "Fraction of Genes Having the Consequence Indicated\n(of Genes Affected by Either of Opposing Consequences)\n(With 95% Confidence Interval)"
consequenceBalance$nTot <- consequenceBalance$nDown + 
  consequenceBalance$nUp
consequenceBalance$Comparison <- paste(consequenceBalance$condition_1, 
                                        "vs", consequenceBalance$condition_2, sep = "\n")
consequenceBalance$feature2 <- gsub(" \\(", "\n(", 
                                     consequenceBalance$feature)
consequenceBalance$feature2 <- factor(consequenceBalance$feature2, 
                                       levels = rev(sort(unique(as.character(consequenceBalance$feature2)))))
g1 <- ggplot(data = consequenceBalance, aes(y = feature2,
                                            x = propOfRelevantEvents, color = Significant)) + 
  geom_errorbarh(aes(xmax = propCiLo, xmin = propCiHi), 
                 height = 0.3) + geom_point(aes(size = nTot)) + 
  facet_wrap(~Comparison) + geom_vline(xintercept = 0.5, 
                                       linetype = "dashed") + labs(x = xText, y = "Consequence of Isoform Switch\n(and the opposing consequence)") + 
  theme_bw() + theme(axis.text.x = element_text(angle = -45, 
                                                hjust = 0, vjust = 1)) + scale_color_manual(name = paste0("FDR < ", 
                                                                                                          0.05), values = c("black", "red"), drop = FALSE) + 
  scale_radius(limits = c(0, max(roundUpToNearsTenOrHundred(consequenceBalance$nTot)))) + 
  guides(color = guide_legend(order = 1), size = guide_legend(order = 2)) + 
  coord_cartesian(xlim = c(0, 1))
g1 <- g1 + labs(size = "Genes")
g1
ggsave(plot = g1,'Result/Functional_consequences_enrichment.pdf',width =10,height = 6)
