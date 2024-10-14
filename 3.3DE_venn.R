library(ggplot2)
library(tidyverse)
library(Seurat)
library(ggVennDiagram)
id_map = read.table('Result/id_map.txt',sep = '\t',header = T)
#DE_venn
## dominant tr 
dominant_transcripts_pct_ad_normal = readRDS('Result/dominant_transcripts_pct_ad_normal.RDS')

dominant_transcripts_pct_ad_normal$lgFC = log2((dominant_transcripts_pct_ad_normal$ad_pct+0.000001)/(dominant_transcripts_pct_ad_normal$nor_pct+0.000001))
dominant_transcripts_pct_ad_normal = na.omit(dominant_transcripts_pct_ad_normal)
dominant_transcripts_pct_ad_normal$type = ifelse(dominant_transcripts_pct_ad_normal$lgFC>1,'AD',
                                                 ifelse(dominant_transcripts_pct_ad_normal$lgFC < -1 ,'Normal','notSig'))
table(dominant_transcripts_pct_ad_normal$type)
dominant_transcripts_pct_ad_normal = left_join(dominant_transcripts_pct_ad_normal,unique(id_map[,c('gene_id','gene_name')]))
dominant_tr_change_ge = unique(dominant_transcripts_pct_ad_normal$gene_name[dominant_transcripts_pct_ad_normal$type !='notSig'] )

#DTE
sce.nano = readRDS('Result/AD_Keratinocytes_nanopore.rds')
DefaultAssay(sce.nano) = 'RNA'
Idents(sce.nano) = sce.nano$stage
DTE = FindMarkers(sce.nano,ident.1 = 'AD',ident.2 = 'Normal',logfc.threshold = 0)
DTE$change = ifelse(DTE$p_val_adj <0.05 & DTE$avg_log2FC > log2(1.5),'Up',
                    ifelse(DTE$p_val_adj < 0.05 & DTE$avg_log2FC < -log2(1.5),'Down','Stable'))
table(DTE$change)
DTE$transcript_id = rownames(DTE)
DTE = left_join(DTE,id_map)
DTE_gene = DTE$gene_name[DTE$change != 'Stable']

#DEG
sce.ill = readRDS('Result/AD_Keratinocytes.rds')
DefaultAssay(sce.ill) = 'RNA'
Idents(sce.ill) = sce.ill$stage
DGE = FindMarkers(sce.ill,ident.1 = 'AD',ident.2 = 'Normal',logfc.threshold = 0)
DGE$change = ifelse(DGE$p_val_adj <0.05 & DGE$avg_log2FC > log2(1.5),'Up',
                    ifelse(DGE$p_val_adj < 0.05 & DGE$avg_log2FC < -log2(1.5),'Down','Stable'))
table(DGE$change)
DGE_gene = rownames(DGE)[DGE$change != 'Stable']

#gene pct diff
# save.image('Result/DE_venn.rda')
# load('Result/DE_venn.rda')
gene_pct = read_rds('Result/gene_pct_ad_normal.RDS')
gene_pct$lgFC = log2((gene_pct$ad_pct+0.000001)/(gene_pct$nor_pct+0.000001))
gene_pct = na.omit(gene_pct)
gene_pct$type = ifelse(gene_pct$lgFC>1,'AD',
                                                 ifelse(gene_pct$lgFC < -1 ,'Normal','notSig'))
table(gene_pct$type)
ge_pct_change_ge = gene_pct$gene_name[gene_pct$type != 'notSig']


####PLOT

gene_list <- list(
  "DGU (Differential Gene Usage)" = ge_pct_change_ge,
  "DGE (Differential Gene Expression)" = DGE_gene,
  "DTE (Differential Transcript Expression)" = DTE_gene,
  "DDTU (Differential Dominant Transcript Usage)" = dominant_tr_change_ge
)


ggVennDiagram(gene_list) +
  scale_fill_gradient(low = "white", high = "steelblue")

ggsave('Result/DE_venn.pdf',width = 8,height = 6)
