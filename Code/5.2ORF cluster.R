library(Seurat)
library(harmony)
library(readr)
library(data.table)
library(clustree)
library(tidyverse)
sce.ill = read_rds('AD_Keratinocytes.rds')
sce.nano = read_rds('AD_Keratinocytes_nanopore_gene_exp_seurat.rds')
common_cell = intersect(colnames(sce.ill),colnames(sce.nano))
orf_exp = fread('orf_mer_exp.csv',sep = ',',header = T,data.table = F)
rownames(orf_exp) = orf_exp$orf_name
orf_exp_sub = orf_exp[,common_cell]
sce.orf = CreateSeuratObject(orf_exp_sub,meta.data = sce.ill@meta.data)
sce.orf <- NormalizeData(sce.orf)
sce.orf <- FindVariableFeatures(sce.orf)
sce.orf <- ScaleData(sce.orf, verbose=FALSE) # normally do not need regress out the sources of heterogeneity from the data
sce.orf <- RunPCA(sce.orf, npcs=30,verbose=FALSE)
sce.orf <- RunHarmony(sce.orf,'orig.ident')
sce.orf <- FindNeighbors(sce.orf, reduction="harmony",dims = 1:30)
sce.orf <- FindClusters(sce.orf, resolution = c(seq(.1,1,.1)))
for (i in seq(.1,1,.1)) {
  sce.orf@meta.data[,paste0('RNA_snn_res.',i)][sce.orf@meta.data[,paste0('RNA_snn_res.',i)] == max(levels(sce.orf@meta.data[,paste0('RNA_snn_res.',i)]))] = levels(sce.orf@meta.data[,paste0('RNA_snn_res.',i)])[length(levels(sce.orf@meta.data[,paste0('RNA_snn_res.',i)]))-1]
}
pdf('orf_clustree.pdf',width = 8,height = 6)
clustree(sce.orf@meta.data, prefix = "RNA_snn_res.")
dev.off()
sce.orf <- RunUMAP(sce.orf,reduction="harmony",dims=1:30)
Idents(sce.orf) = sce.orf$RNA_snn_res.0.4
p1 = DimPlot(sce.orf,label = T)
p1
saveRDS(sce.orf,'AD_Keratinocytes_nanopore_orf_exp_seurat.rds')
markers = FindAllMarkers(sce.orf,only.pos = T)
marker_sig <- markers %>%
  dplyr::filter(abs(avg_log2FC)>log2(1.5), p_val_adj<0.05, pct.1 > 0.2)
top_markers <- marker_sig %>%
  group_by(cluster) %>%
  top_n(5,wt = avg_log2FC)
marker_sig_6 =  markers %>%
  dplyr::filter(cluster == 6,abs(avg_log2FC)>log2(1.5), p_val_adj<0.05, pct.1 > 0.08)
m_0 = c('ORF-119462','ORF-32916','ORF-25217','ORF-109851','ORF-14435')
m_1 = top_markers$gene[top_markers$cluster=='1']
m_2 = top_markers$gene[top_markers$cluster=='2']
m_3 = top_markers$gene[top_markers$cluster=='3']
m_4 = top_markers$gene[top_markers$cluster=='4']
m_5 = top_markers$gene[top_markers$cluster=='5']
m_6 = marker_sig_6$gene[1:5]

top_markers_final = data.frame(cluster = rep(0:6,each = 5),orf = c(m_0,m_1,m_2,m_3,m_4,m_5,m_6))
top_markers_final$orf = gsub('-','_',top_markers_final$orf)
top_markers_final = left_join(top_markers_final,orf_df,by = c('orf'='orf_name'))
orf_df = fread('ORF_id_map.txt',header = T,sep = '\t')
orf_df_sub = orf_df[orf_df$orf_name %in% gsub('-','_',top_markers_final$orf),]
write.table(orf_df_sub,'AD_Keratinocytes_orf_marker_dotplot_seq.txt',sep = '\t',row.names = F,quote = F)
p <- DotPlot(object = sce.orf, features = unique(top_markers_final$orf)) +
  scale_color_gradient2(low = "blue", mid = "#D9D9D9", high = "red") +
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = .5))
ggsave(plot = p,filename = 'AD_Keratinocytes_orf_marker_dotplot.pdf',width = 8,height = 4)


library(ggsankey)
cell_num_df = as.data.frame(table(Idents(sce.orf)))
colnames(cell_num_df) = c('cluster','num')
cell_num_df$cluster_with_num = paste0(cell_num_df$cluster,' (',cell_num_df$num,')')
meta = sce.orf@meta.data
meta$cell_id = rownames(meta)
meta = left_join(meta,cell_num_df,by =c('RNA_snn_res.0.4'='cluster'))
identical(rownames(sce.orf@meta.data),meta$cell_id)
sce.orf$cluster_with_num = meta$cluster_with_num
pdf('/data/haowu/sc_long/subclass/UMAP_orf_cluster.pdf',width = 8,height = 6)
DimPlot(sce.orf,group.by = 'cluster_with_num')
dev.off()


#meta_ill = sce.ill@meta.data[common_cell,]
common_cell = intersect(colnames(sce.nano),colnames(sce.orf))
meta_nano = sce.nano@meta.data[common_cell,]
meta_orf = sce.orf@meta.data[common_cell,]
cp_dt = data.frame(nano = meta_nano$RNA_snn_res.0.3,
                   orf = meta_orf$RNA_snn_res.0.4)
all_ctypes <- factor(c(unique(cp_dt$nano),unique(cp_dt$orf)))

ordered_ctypes <-levels(all_ctypes)

ctype_map <- tibble(name = ordered_ctypes, id = 1:length(unique(all_ctypes)))
# 
# ashok_data_before <- inner_join(tibble(name = cp_dt$ill), ctype_map)
ashok_data_before <- inner_join(tibble(name = cp_dt$nano), ctype_map)
ashok_data_after <- inner_join(tibble(name = cp_dt$orf), ctype_map)
alluvial_data <- tibble(before = ashok_data_before$id,
                        after = ashok_data_after$id)

df <- alluvial_data %>%
  make_long(before,after)

node_name_tbl <- tibble(node = ctype_map$id, name = ctype_map$name)
df <- inner_join(df, node_name_tbl)
df <- df %>% arrange(node, next_node)

pdf('nano_ge_orf_Sankey.pdf',width = 6,height = 10)
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
dev.off()

# cp_dt_sub = cp_dt[cp_dt$nano == 'Keratinocytes',]
# data = as.data.frame(table(cp_dt_sub$ill))
# colnames(data) = c('Celltype','Number')
# data <- data %>%
#   mutate(Percentage = Number / sum(Number) * 100,
#          legend_label = paste0(Celltype, " (", round(Percentage, 1), "%)"))
# 
# # Plot
# p2 = ggplot(data, aes(x = "", y = Percentage, fill = legend_label)) +
#   geom_bar(stat = "identity", width = 1) +
#   coord_polar("y", start = 0) +
#   theme_void() +
#   theme(legend.position = "right") +
#   guides(fill = guide_legend(title = "Cell Type")) +
#   scale_fill_brewer(palette = "Paired")
# print(p2)
# dev.off()
# 
pdf('nano_ge_orf_pie_in_nano.pdf',width = 8,height = 6)
for (ct in unique(cp_dt$nano)) {
  cp_dt_sub = cp_dt[cp_dt$nano == ct,]
  data = as.data.frame(table(cp_dt_sub$orf))
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


pdf('nano_ge_orf_pie_in_orf.pdf',width = 8,height = 6)
for (ct in unique(cp_dt$orf)) {
  cp_dt_sub = cp_dt[cp_dt$orf == ct,]
  data = as.data.frame(table(cp_dt_sub$nano))
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


