library(UCell)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(readr)
library(harmony)
library(clustree)
#https://doi.org/10.1016/j.jaci.2023.04.019
Natalia_KC1 = c('KRT15','COL17A1','S100A2','KRT5','DST')
Natalia_KC2 = c('KRT2','KRT10','KRT16','KRTDAP','KRT1')
Natalia_KC3 = c('AQP3','LGALS7','LGALS7B','S100A2','CCL27')
Natalia_KC4 = c('CALML3','DEFB1','MIR205HG','GPX2','KRT17')
Natalia_KCpro = c('KRT5','STMN1','S100A2','KRT14','PTTG1')
#https://doi.org/10.1016/j.jaci.2020.01.042
Helen_Basal = c('KRT15','GAS5','IL1R2','DST','SERPINF1')
Helen_Channel_ATPase = c('DEFB1','TFCP2L1','C19orf33','SEMA3C','ATP1B3')
Helen_Channel_Gap = c('BNC2','LGR5','CREB5','LHX2','COMP')
Helen_IRS_Sebaceous = c('FASN','TF','KRT79','ALOX15B','ACSBG1')
Helen_ORS = c('KRT6B','CST6','KRT17','SBSN','KRT6A')
Helen_Proliferating = c('UBE2C','TOP2A','CENPF','MKI67','ASPM')
Helen_Suprabasal = c('KRTDAP','SERPINB4','KRT16','KRT1','DMKN')
#https://doi.org/10.1016/j.jaci.2020.03.041
Thomas_KC1 = c('LGALS7','NUPR1','KLF6','MT1X','CSTA')
Thomas_KC2 = c('RPS18','TMSB4X','KLF6','AVPI1','CSTA')
Thomas_KC3 = c('LGALS7','TMSB4X','NUPR1','TPPP3','CSTA')
Thomas_KC4 = c('LGALS7','TPPP3','CSTA','HSPB1','LYPD3')
Thomas_KC5 = c('LGALS7','LYPD3','CSTA','TPPP3','KLF6')
Thomas_KC6 = c('CRABP2','KLF6','CSTA','LGALS7','KRT1')
Thomas_KC7 = c('RPS18','MT1X','NUPR1','TPPP3','LGALS7')
#https://doi.org/10.1111/jdv.19256
Jie_Basal = c('KRT15','COL17A1','POSTN','CXCL14','IGFBP3')
Jie_Suprabasal = c('KRT1','KRT10','KRTDAP','DMKN','KRT2')
Jie_Proliferating = c('STMN1','TUBA1B','PTTG1','HMGB2','UBE2C')
Jie_Late_differentiated = c('SPRR2E','SPRR2G','LCE3D','CNFN','FLG')

signature = list(Natalia_KC1,Natalia_KC2,Natalia_KC3,Natalia_KC4,Natalia_KCpro,Helen_Basal,Helen_Channel_ATPase, Helen_Channel_Gap, Helen_IRS_Sebaceous, Helen_ORS, Helen_Proliferating, Helen_Suprabasal,Thomas_KC1,Thomas_KC2,Thomas_KC3,Thomas_KC4,Thomas_KC5,Thomas_KC6,Thomas_KC7,Jie_Basal,Jie_Suprabasal,Jie_Proliferating,Jie_Late_differentiated)
names(signature) = c('Natalia_KC1','Natalia_KC2','Natalia_KC3','Natalia_KC4','Natalia_KCpro','Helen_Basal','Helen_Channel_ATPase', 'Helen_Channel_Gap', 'Helen_IRS_Sebaceous', 'Helen_ORS', 'Helen_Proliferating', 'Helen_Suprabasal','Thomas_KC1','Thomas_KC2','Thomas_KC3','Thomas_KC4','Thomas_KC5','Thomas_KC6','Thomas_KC7','Jie_Basal','Jie_Suprabasal','Jie_Proliferating','Jie_Late_differentiated')
sce.ill = read_rds('Result/AD_Keratinocytes.rds')
DefaultAssay(sce.ill) = 'RNA'
markers = FindAllMarkers(sce.ill,only.pos = T)
marker_sig <- markers %>% 
  filter(abs(avg_log2FC)>log2(1.5),p_val_adj<0.05)
marker_sig = marker_sig[marker_sig$gene %in% rownames(sce.nano.ge),]
top_markers <- marker_sig %>% 
  group_by(cluster) %>% 
  top_n(5,wt = avg_log2FC)
top_markers$cluster = paste0('paired_KC',top_markers$cluster)
paired_sig_list = split(top_markers$gene,top_markers$cluster)

signature = c(signature,paired_sig_list)
save(signature,file = '/data/haowu/sc_long//signature_cor/signature.Rdata')
load('/data/haowu/sc_long//signature_cor/signature.Rdata')
sce.nano.ge_exp = read_rds('Result/AD_Keratinocytes_nanopore_gene_exp.rds')
sce.nano.ge_exp = as.data.frame(sce.nano.ge_exp)
rownames(sce.nano.ge_exp) = sce.nano.ge_exp$gene_name
sce.nano.ge_exp$gene_name = NULL
common_cells = intersect(colnames(sce.ill),colnames(sce.nano.ge_exp))
sce.nano.ge_exp = sce.nano.ge_exp[,common_cells]
sce.nano.ge = CreateSeuratObject(sce.nano.ge_exp,meta.data = sce.ill@meta.data)
sce.nano.ge <- NormalizeData(sce.nano.ge)
sce.nano.ge <- FindVariableFeatures(sce.nano.ge)
sce.nano.ge <- ScaleData(sce.nano.ge, verbose=FALSE) # normally do not need regress out the sources of heterogeneity from the data
sce.nano.ge <- RunPCA(sce.nano.ge, npcs=30,verbose=FALSE)
sce.nano.ge <- RunHarmony(sce.nano.ge,'orig.ident')
sce.nano.ge <- FindNeighbors(sce.nano.ge, reduction="harmony",dims = 1:30)
sce.nano.ge <- FindClusters(sce.nano.ge, resolution = seq(0.1,0.8,0.05))
pdf("Result/AD_Keratinocytes_nanopore_gene_exp_marker_check_dotplot.pdf", height = 8, width = 14)
for (i in seq(0.1,0.8,0.05)) {
  Idents(sce.nano.ge) = sce.nano.ge@meta.data[[paste0("RNA_snn_res.",i)]]
  markers = FindAllMarkers(sce.nano.ge,only.pos = T)
  marker_sig <- markers %>% 
    dplyr::filter(abs(avg_log2FC)>log2(1.5),p_val_adj < 0.05,pct.1>0.2)
  marker_sig =marker_sig[!grepl('HLA',marker_sig$gene),]
  top_markers <- marker_sig %>% 
    group_by(cluster) %>% 
    top_n(5,wt = avg_log2FC)
  p <- DotPlot(object = sce.nano.ge, features = unique(top_markers$gene)) +
    scale_color_gradient2(low = "blue", mid = "#D9D9D9", high = "red") +
    theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = .5))
  
  print(p)
  
}
dev.off()

pdf('Result/AD_Keratinocytes_nanopore_gene_exp_seurat_clustree.pdf',width = 10,height = 8)
clustree(sce.nano.ge@meta.data, prefix = "RNA_snn_res.")
dev.off()
sce.nano.ge <- RunUMAP(sce.nano.ge,reduction="harmony",dims=1:30)
Idents(sce.nano.ge) = sce.nano.ge$RNA_snn_res.0.3
saveRDS(sce.nano.ge,'Result/AD_Keratinocytes_nanopore_gene_exp_seurat.rds')
DimPlot(sce.nano.ge,label = T)
ggsave('Result/AD_Keratinocytes_nanopore_gene_exp_UMAP_v2_res0.3.pdf',width = 8,height = 6)
DefaultAssay(sce.nano.ge)
markers = FindAllMarkers(sce.nano.ge,only.pos = T)
marker_sig <- markers %>% 
  dplyr::filter(abs(avg_log2FC)>log2(1.5),p_val_adj < 0.05,pct.1>0.2)
marker_sig =marker_sig[!grepl('HLA',marker_sig$gene),]
top_markers <- marker_sig %>% 
  group_by(cluster) %>% 
  top_n(5,wt = avg_log2FC)

# ####
p <- DotPlot(object = sce.nano.ge, features = unique(top_markers$gene)) +
  scale_color_gradient2(low = "blue", mid = "#D9D9D9", high = "red") +
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = .5))
pdf("Result/AD_Keratinocytes_nanopore_gene_exp_marker_dotplot.pdf", height = 8, width = 14)
print(p)
dev.off()


signature.matrix.ker <- (data.frame(ScoreSignatures_UCell(as.data.frame(sce.nano.ge@assays$RNA$data), features = signature,ncores= 20)))
# signature.matrix.ker <- (data.frame(ScoreSignatures_UCell(as.data.frame(sce.nano.ge@assays$RNA$data), features = paired_sig_list,ncores= 20)))
save(signature.matrix.ker,file = '/data/haowu/sc_long//signature_cor/signature_exp.Rdata')
sce.nano.ge@meta.data = cbind(sce.nano.ge@meta.data, signature.matrix.ker)
colnames(sce.nano.ge@meta.data )
sce.nano.ge$cellIdent_sub = factor(paste0('Ker_',sce.nano.ge$RNA_snn_res.0.3),levels = paste0('Ker_',0:7))
Idents(sce.nano.ge) = sce.nano.ge$cellIdent_sub
DotPlot(sce.nano.ge, features = colnames(signature.matrix.ker)) +
  ggplot2::coord_flip() +
  RotatedAxis() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0)
ggsave('/data/haowu/sc_long/signature_cor/Comparison_published_Dotplot_v3_res0.3.pdf',width = 8,height = 10)

#ill signature
markers_nano_ge = FindAllMarkers(sce.nano.ge,only.pos = T)
marker_sig_nano_ge <- markers_nano_ge %>% 
  filter(abs(avg_log2FC)>log2(1.5),p_val_adj<0.05)
marker_sig_nano_ge = marker_sig_nano_ge[marker_sig_nano_ge$gene %in% rownames(sce.ill),]
top_markers <- marker_sig_nano_ge %>% 
  group_by(cluster) %>% 
  top_n(5,wt = avg_log2FC)
top_markers$cluster = paste0('nano_',top_markers$cluster)
paired_sig_list = split(top_markers$gene,top_markers$cluster)
signature = c(signature,paired_sig_list)
sce.ill = read_rds('Result/AD_Keratinocytes.rds')
DefaultAssay(sce.ill) = 'RNA'

signature.matrix.ker_ill <- (data.frame(ScoreSignatures_UCell(as.data.frame(sce.ill@assays$RNA$data), features = signature,ncores= 20)))  
sce.ill@meta.data = cbind(sce.ill@meta.data, signature.matrix.ker_ill)
colnames(sce.ill@meta.data)
sce.ill$cellIdent_sub = factor(paste0('ILL_Ker_',sce.ill$seurat_clusters),levels = paste0('ILL_Ker_',0:5))
Idents(sce.ill) = sce.ill$cellIdent_sub
DotPlot(sce.ill, features = colnames(signature.matrix.ker_ill)) +
  ggplot2::coord_flip() +
  RotatedAxis() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0)
ggsave('/data/haowu/sc_long/signature_cor/Comparison_published_Dotplot_ill.pdf',width = 8,height = 10)
