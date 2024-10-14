library(tidyverse)
library(Seurat)
library(ggsankey)
library(readr)
library(harmony)
library(clustree)
library(SCINA)
sce.ill = read_rds('Result/AD_Keratinocytes.rds')
sce.nano = read_rds('Result/AD_Keratinocytes_nanopore.rds')  
DefaultAssay(sce.ill) = 'RNA'
DefaultAssay(sce.nano) = 'RNA'
ill_exp = sce.ill@assays$RNA@counts[,colnames(sce.nano)]
sce.ill = CreateSeuratObject(ill_exp,meta.data = sce.ill@meta.data)
sce.ill <- NormalizeData(sce.ill)
sce.ill <- FindVariableFeatures(sce.ill)
sce.ill <- ScaleData(sce.ill, verbose=FALSE) # normally do not need regress out the sources of heterogeneity from the data
sce.ill <- RunPCA(sce.ill, npcs=30,verbose=FALSE)
sce.ill <- RunHarmony(sce.ill,'orig.ident')
sce.ill <- FindNeighbors(sce.ill, reduction="harmony",dims = 1:30)
sce.ill <- FindClusters(sce.ill, resolution = c(seq(.1,1,.1)))
clustree(sce.ill@meta.data, prefix = "RNA_snn_res.")
sce.ill <- RunUMAP(sce.ill,reduction="harmony",dims=1:30)
Idents(sce.ill) = sce.ill$RNA_snn_res.0.7
p1 = DimPlot(sce.ill,label = T)

#
sce.nano <- NormalizeData(sce.nano)
sce.nano <- FindVariableFeatures(sce.nano)
sce.nano <- ScaleData(sce.nano, verbose=FALSE) # normally do not need regress out the sources of heterogeneity from the data
sce.nano <- RunPCA(sce.nano, npcs=30,verbose=FALSE)
sce.nano <- RunHarmony(sce.nano,'orig.ident')
sce.nano <- FindNeighbors(sce.nano, reduction="harmony",dims = 1:30)
sce.nano <- FindClusters(sce.nano, resolution = c(seq(.1,1,.1)))
clustree(sce.nano@meta.data, prefix = "RNA_snn_res.")
sce.nano <- RunUMAP(sce.nano,reduction="harmony",dims=1:30)
Idents(sce.nano) = sce.nano$RNA_snn_res.0.6
p2 = DimPlot(sce.nano,label = T)

#
sce.nano.ge_exp = read_rds('Result/AD_Keratinocytes_nanopore_gene_exp.rds')
sce.nano.ge_exp = as.data.frame(sce.nano.ge_exp)
rownames(sce.nano.ge_exp) = sce.nano.ge_exp$gene_name
sce.nano.ge_exp$gene_name = NULL
identical(colnames(sce.nano.ge_exp),colnames(sce.ill))
sce.nano.ge = CreateSeuratObject(sce.nano.ge.ge_exp,meta.data = sce.ill@meta.data)
sce.nano.ge <- NormalizeData(sce.nano.ge)
sce.nano.ge <- FindVariableFeatures(sce.nano.ge)
sce.nano.ge <- ScaleData(sce.nano.ge, verbose=FALSE) # normally do not need regress out the sources of heterogeneity from the data
sce.nano.ge <- RunPCA(sce.nano.ge, npcs=30,verbose=FALSE)
sce.nano.ge <- RunHarmony(sce.nano.ge,'orig.ident')
sce.nano.ge <- FindNeighbors(sce.nano.ge, reduction="harmony",dims = 1:30)
sce.nano.ge <- FindClusters(sce.nano.ge, resolution = c(seq(.1,1,.1)))
clustree(sce.nano.ge@meta.data, prefix = "RNA_snn_res.")
sce.nano.ge <- RunUMAP(sce.nano.ge,reduction="harmony",dims=1:30)
Idents(sce.nano.ge) = sce.nano.ge$RNA_snn_res.0.8
p3 = DimPlot(sce.nano.ge,label = T)

cp_dt = data.frame(nano_ge = sce.nano.ge@meta.data$RNA_snn_res.0.4,
                   nano_tr = sce.nano@meta.data$RNA_snn_res.0.4)
all_ctypes <- c(unique(cp_dt$nano_tr), unique(cp_dt$nano_ge))

ordered_ctypes <- levels(all_ctypes)[order(as.numeric(levels(all_ctypes)))]

ctype_map <- tibble(name = ordered_ctypes, id = 1:length(unique(all_ctypes)))

ashok_data_before <- inner_join(tibble(name = cp_dt$nano_ge), ctype_map)
ashok_data_after <- inner_join(tibble(name = cp_dt$nano_tr), ctype_map)
alluvial_data <- tibble(before = ashok_data_before$id,
                        after = ashok_data_after$id)

df <- alluvial_data %>%
  make_long(before, after)

node_name_tbl <- tibble(node = ctype_map$id, name = ctype_map$name)
df <- inner_join(df, node_name_tbl)
df <- df %>% arrange(node, next_node)
ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = name)) +
  geom_alluvial(flow.alpha = .7, width = 0.25, space = 300) +
  geom_alluvial_label(size = 5, color = "black", fill = "white", space = 300) +
  scale_fill_viridis_d() +
  theme_alluvial() +
  labs(x = NULL, y = NULL) +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
save.image('Result/recluster.rdata')



#####all cell type
sce.ill = read_rds('Resultintegration_Illumina.rds')
sce.nano.ge_exp = data.table::fread('ResultAD_nanopore_data_ge_exp_count.txt',sep = '\t',data.table = F)
common_cell = intersect(colnames(sce.ill),colnames(sce.nano.ge_exp))
rownames(sce.nano.ge_exp) = sce.nano.ge_exp$gene_name
sce.nano.ge_exp$gene_name = NULL
sce.nano.ge_exp = sce.nano.ge_exp[,common_cell]
meta = sce.ill@meta.data[common_cell,]
sce.nano.ge = CreateSeuratObject(sce.nano.ge_exp,meta.data = meta)
sce.nano.ge <- NormalizeData(sce.nano.ge)
sce.nano.ge <- FindVariableFeatures(sce.nano.ge)
sce.nano.ge <- ScaleData(sce.nano.ge, verbose=FALSE) # normally do not need regress out the sources of heterogeneity from the data
sce.nano.ge <- RunPCA(sce.nano.ge, npcs=30,verbose=FALSE)
sce.nano.ge <- RunHarmony(sce.nano.ge,'orig.ident')
sce.nano.ge <- FindNeighbors(sce.nano.ge, reduction="harmony",dims = 1:30)
sce.nano.ge <- RunUMAP(sce.nano.ge,reduction="harmony",dims=1:30)
sce.nano.ge <- FindClusters(sce.nano.ge, resolution = 1)
DimPlot(sce.nano.ge,label = T)
DefaultAssay(sce.nano.ge) <- "RNA"

all_markers <- FindAllMarkers(sce.nano.ge, only.pos = TRUE)
library(SingleR)
refdata <- get(load("/data/haowu/reference/singleR/ref_Human_all.RData"))
testdata <- GetAssayData(sce.nano.ge, slot="data")
clusters <- sce.nano.ge@meta.data$seurat_clusters
cellpred <- SingleR(test = testdata, ref = refdata, 
                    labels = refdata$label.fine, clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")
table(cellpred$labels)
cell.num=as.matrix(table(sce.nano.ge@meta.data$seurat_clusters))
celltype <- data.frame(ClusterID=rownames(cellpred),
                       celltype=cellpred$labels, stringsAsFactors = F)
celltype$celltype.main=gsub(":.*","",celltype$celltype)

marker_collected <- read.table('Resultcell_marker_collected.txt',sep = "\t",header = T,fill = T)
# markers <- c("KRT15","KRT1","TOP2","AKRT14","KRT5","AQP5","DCD","SELE","MYL9","IL22","IL13","CCL5","NKG7","GZMB","DUT","TYMS","CTLA4","TNFRSF18","TIGIT","TPSAB1","TPSB2","IL1B","IRF8","LAMP3","CCL22","CCL17","CD207","CD1A","F13A1","CCL18","CCL13","MLANA","PMEL","KIT","SOX10","TYR","MITF","COL1A1")
markers <- strsplit(marker_collected$markers,',')
names(markers) = marker_collected$cell
p1 = DotPlot(object = sce.nano.ge, features = markers) +
  scale_color_gradient2(low = "blue", mid = "#D9D9D9", high = "red") +
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = .5))
pdf( "ResultAD_key_cluster_markers.pdf", height = 8, width = 30)
print(p1)
dev.off()
p2 = DimPlot(sce.nano.ge,label = T)
p1/p2
cluster_type <- c("Keratinocytes","Keratinocytes","Keratinocytes","Keratinocytes","Keratinocytes","Keratinocytes","Keratinocytes","Keratinocytes","Keratinocytes","Keratinocytes","Natural_killer_cells","T_cells","Keratinocytes","Smooth_muscle_cells","Fibroblasts","Melanocytes","Sweat_glands_cells")
manual_annotation = data.frame(ClusterID = levels(sce.nano.ge$seurat_clusters),
                               celltype.main = cluster_type)
sce.nano.ge@meta.data$manual_annotation = "Unknown"
for (i in 1:nrow(manual_annotation)) {
  sce.nano.ge@meta.data[which(sce.nano.ge@meta.data$seurat_clusters == manual_annotation$ClusterID[i]),'manual_annotation'] <- manual_annotation$celltype.main[i]
}
table(sce.nano.ge@meta.data$manual_annotation)
DimPlot(sce.nano.ge,group.by = 'manual_annotation',label = T)
ggsave('ResultUMAP_nano_ge_manual_annotation.pdf',width = 8,height = 6)
# tmp = sce.ill@meta.data[rownames(sce.nano.ge@meta.data)[sce.nano.ge$seurat_clusters == 14],]
# tmp = sce.nano.ge@meta.data[rownames(sce.ill@meta.data)[sce.ill$manual_annotation == 'Natural_killer_cells'],]
# table(tmp$seurat_clusters)
library(ggsankey)
meta_ill = sce.ill@meta.data[common_cell,]
meta_nano = sce.nano.ge@meta.data[common_cell,]
cp_dt = data.frame(ill = meta_ill$manual_annotation,
                   nano = meta_nano$manual_annotation)
all_ctypes <- factor(c(unique(cp_dt$ill), unique(cp_dt$nano)))

ordered_ctypes <-levels(all_ctypes)

ctype_map <- tibble(name = ordered_ctypes, id = 1:length(unique(all_ctypes)))

ashok_data_before <- inner_join(tibble(name = cp_dt$nano), ctype_map)
ashok_data_after <- inner_join(tibble(name = cp_dt$ill), ctype_map)
alluvial_data <- tibble(before = ashok_data_before$id,
                        after = ashok_data_after$id)

df <- alluvial_data %>%
  make_long(before, after)

node_name_tbl <- tibble(node = ctype_map$id, name = ctype_map$name)
df <- inner_join(df, node_name_tbl)
df <- df %>% arrange(node, next_node)

pdf('Resultill_nano_ge_Sankey.pdf',width = 6,height = 10)
p1 = ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = name)) +
  geom_alluvial(flow.alpha = .7, width = 0.25, space = 300) +
  geom_alluvial_label(size = 3, color = "black", fill = "white", space = 300) +
  scale_fill_brewer(palette = "Paired") +
  theme_alluvial() +
  labs(x = NULL, y = NULL) +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
print(p1)
cp_dt_sub = cp_dt[cp_dt$nano == 'Keratinocytes',]
data = as.data.frame(table(cp_dt_sub$ill))
colnames(data) = c('Celltype','Number')
data <- data %>%
  mutate(Percentage = Number / sum(Number) * 100,
         legend_label = paste0(Celltype, " (", round(Percentage, 1), "%)"))

# Plot
p2 = ggplot(data, aes(x = "", y = Percentage, fill = legend_label)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  theme_void() +
  theme(legend.position = "right") +
  guides(fill = guide_legend(title = "Cell Type")) +
  scale_fill_brewer(palette = "Paired")
print(p2)
dev.off()

pdf('Resultill_nano_ge_pie.pdf',width = 8,height = 6)
for (ct in unique(cp_dt$nano)) {
  cp_dt_sub = cp_dt[cp_dt$nano == ct,]
  data = as.data.frame(table(cp_dt_sub$ill))
  colnames(data) = c('Celltype','Number')
  data <- data %>%
    mutate(Percentage = Number / sum(Number) * 100,
           legend_label = paste0(Celltype, " (", round(Percentage, 1), "%)"))
  
  # Plot
  p = ggplot(data, aes(x = "", y = Percentage, fill = legend_label)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    theme_void() +
    theme(legend.position = "right") +
    guides(fill = guide_legend(title = "Cell Type")) +
    scale_fill_brewer(palette = "Paired")+labs(title = ct)
  print(p)
}
dev.off()
