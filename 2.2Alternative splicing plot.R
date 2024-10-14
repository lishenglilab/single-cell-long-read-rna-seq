# 加载必要的库
library(data.table)
library(dplyr)
library(ggplot2)
library(ggsci)
library(ggthemes)
library(pheatmap)
library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)


ASE <- fread("Result/merged.0.05filter.psi", sep = "\t", data.table = F)
sce_sub <- readRDS("Result/AD_Keratinocytes.rds")


ASE <- ASE[, c('ASE', intersect(colnames(ASE), colnames(sce_sub)))]
info <- data.frame(cell_id = colnames(ASE)[-1])
info <- left_join(info, sce_sub@meta.data)


res <- data.frame()
res_ls <- vector('list')

for (i in unique(info$seurat_clusters)) {
  cells <- info[info$seurat_clusters == i,]$cell_id
  ASE_sub <- ASE[, c("ASE", cells)]
  rownames(ASE_sub) <- ASE_sub$ASE
  ASE_sub$ASE <- NULL
  count <- apply(ASE_sub, 1, function(x) sum(x > 0.1))
  select_ase <- names(count)[count >= 0.05 * (ncol(ASE) - 1)]
  res <- rbind(res, data.frame(ASE = select_ase, k_sub = i))
  res_ls[[i]] <- select_ase
}
res$type <- gsub(".*;([^:]+):.*", "\\1", res$ASE)

plt_data <- res %>% group_by(k_sub, type) %>% summarise(count = n())
ggplot(plt_data, aes(k_sub, type, fill = count)) +
  geom_tile() +
  scale_fill_binned_sequential(palette = "Purples 3", begin = 0.5, end = 1) +
  geom_text(aes(label = count)) +
  theme_base()
ggsave('Result/AS_ksub_heatmap.pdf', width = 8, height = 6)


pdf('Result/AS_ksub_upset.pdf', width = 8, height = 6)
upset(fromList(res_ls))
dev.off()

psi <- fread("Result/merged.0.05filter.psi", sep = "\t", data.table = F)
sce_sub <- read_rds("Result/AD_Keratinocytes_nanopore.rds")
col_data <- data.frame(cell_id = colnames(psi)[-1])
col_data$sample <- gsub('.*_', '', col_data$cell_id)
col_data <- left_join(col_data, sce_sub@meta.data[, c('cell_id', 'k_sub_ill')])
col_data$group <- ifelse(str_detect(col_data$sample, 'N'), 'Normal', 'AD')

res <- data.frame()
for (i in 1:nrow(psi)) {
  ase <- psi$ASE[i]
  nor_psi <- psi[i, col_data[col_data$group == 'Normal', ]$cell_id]
  ad_psi <- psi[i, col_data[col_data$group == 'AD', ]$cell_id]
  p <- wilcox.test(as.numeric(nor_psi), as.numeric(ad_psi))$p.value
  fc <- mean(as.numeric(ad_psi)) / mean(as.numeric(nor_psi))
  res <- rbind(res, data.frame(ASE = ase, mean_normal = mean(as.numeric(nor_psi)), mean_ad = mean(as.numeric(ad_psi)), foldchange = fc, p.value = p))
}

res$fdr <- p.adjust(res$p.value, method = 'fdr')
res$change <- ifelse(res$fdr < 0.05 & res$foldchange > 0, 'Up', ifelse(res$fdr < 0.05 & res$foldchange < 0, 'Down', 'Stable'))
write.table(res, 'Result/diff_ase_Keratinocytes.txt', sep = '\t', row.names = F, quote = F)

ase_diff <- fread('Result/diff_ase_Keratinocytes.txt', sep = '\t', header = T)
data <- as.data.frame(sce_sub@reductions$umap@cell.embeddings)
data <- left_join(data, sce_sub@meta.data)
ase_top <- ase_diff %>% filter(change != 'Stable') %>% group_by(change) %>% top_n(5, abs(logFC))

pdf('Result/UMAP_diff_ASE.pdf', width = 12, height = 6)
for (as in unique(ase_top$ASE)) {
  ase_sub <- ase[ase$ASE %in% as,]
  ase_sub <- pivot_longer(ase_sub, cols = colnames(ase_sub)[colnames(ase_sub) != 'ASE'], names_to = 'cell_id', values_to = 'psi')
  tmp <- left_join(data, ase_sub)
  tmp <- na.omit(tmp)
  p <- ggplot(tmp, aes(UMAP_1, UMAP_2, color = psi)) +
    geom_point() +
    scale_color_material() +
    facet_grid(. ~ stage) +
    labs(title = paste0(as, ' (', unique(ase_top$change[ase_top$ASE == as]), ')')) +
    theme_bw()
  print(p)
}
dev.off()

volcano(de = ase_diff, gene_colname = 'ASE', fc_cutoff = 1.5, save_path = 'Result/diff_volcano.pdf')
