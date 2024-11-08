library(ggplot2)
library(dplyr)
library(formattable)
library(stringr)
library(data.table)
library(ggpubr)
library(tidyr)
library(ggsci)
library(ggtranscript)
library(rtracklayer)
sce.nano <- readRDS("integration_nanopore.rds")
sce.ill <- readRDS("integration_Illumina.rds")
DefaultAssay(sce.nano) <- "RNA"
transcript_assembly_classification <- read.table("new_AD_classification.txt",sep = "\t",header = T)
transcript_assembly_classification <- transcript_assembly_classification[transcript_assembly_classification$isoform %in% rownames(sce.nano),]
transcript_assembly_classification$structural_category [!(transcript_assembly_classification$structural_category %in% c("full-splice_match","incomplete-splice_match","novel_in_catalog","novel_not_in_catalog"))] <- "Other"
length(intersect(transcript_assembly_classification$isoform,rownames(sce.nano)))
length(setdiff(transcript_assembly_classification$isoform,rownames(sce.nano)))
id_map = fread('~/project/single_cell_full_length/Result/FLAMES/output/merge/id_map.txt',header = T,sep = '\t')
id_map = id_map[id_map$transcript_id %in% rownames(sce.nano),]
two_tr_ges = names(table(id_map$gene_id)[table(id_map$gene_id) == 2])
two_tr_ges = id_map[id_map$gene_id %in% two_tr_ges,]
cells_stat_df = data.frame()
for (ge in unique(two_tr_ges$gene_id)) {
  trs = two_tr_ges[two_tr_ges$gene_id == ge,]$transcript_id
  exp = as.data.frame(sce.nano@assays$RNA@counts[trs,])
  cells_stat = list()
  for (tr in trs) {
    cells = colnames(exp)[exp[tr,]>0]
    cells_stat[[tr]] = cells
  }
  common_len = length(intersect(cells_stat[[1]],cells_stat[[2]]))
  cells_stat_sub = data.frame(gene_id = ge,
                              all_cell_num = length(unique(unlist(cells_stat))),
                              tr_1_id = trs[1],
                              tr_1_num = length(cells_stat[[1]])-common_len,
                              tr_1_ration = (length(cells_stat[[1]])-common_len)/length(unique(unlist(cells_stat))),
                              tr_2_id = trs[2],
                              tr_2_num = length(cells_stat[[2]])-common_len,
                              tr_2_ration = (length(cells_stat[[2]])-common_len)/length(unique(unlist(cells_stat))),
                              common_num = common_len,
                              common_ration = common_len/length(unique(unlist(cells_stat))))
  cells_stat_df = rbind(cells_stat_df,cells_stat_sub)
}
write.table(cells_stat_df,'two_tr_genes.txt',sep = '\t',row.names = F,quote = F)
gtf <- rtracklayer::import('isoform_annotated.merged.sorted.gtf')
gtf = as.data.frame(gtf)
gtf$transcript_id = gsub('_','-',gtf$transcript_id)
cells_stat_df = read.table('two_tr_genes.txt',sep = '\t',header = T)
cells_stat_df_select = cells_stat_df[cells_stat_df$tr_1_ration>0.35&cells_stat_df$tr_2_ration>0.35 & cells_stat_df$all_cell_num>500,]
cells_stat_df_select = left_join(cells_stat_df_select,unique(id_map[,c('gene_id','gene_type')]))
genes = cells_stat_df_select[cells_stat_df_select$gene_type == 'Protein coding',]$gene_id
genes = c('ENSG00000214717.12_PAR_Y','ENSG00000183971.10','ENSG00000155980.13','ENSG00000279086.1')
cell.embeddings = sce.ill@reductions[["umap"]]@cell.embeddings

pdf('two_example_featureplot.pdf',width = 20,height = 15)
for (ge in genes) {
  trs = id_map[id_map$gene_id %in% ge, ]$transcript_id
  ge_exons = dplyr::filter(gtf, transcript_id %in% trs &type == "exon")
  ge_exons$transcript_id = factor(ge_exons$transcript_id,levels = trs)
  p = ggplot(ge_exons,aes(xstart = start,xend = end,y = transcript_id
  )) +
    geom_range()+
    geom_intron(
      data = to_intron(ge_exons, "transcript_id"),
      aes(strand = strand)
    )+
    theme_bw()
  exp = as.data.frame(sce.nano@assays$RNA@data[trs, ])
  common_cell = intersect(rownames(cell.embeddings),colnames(exp))
  exp = exp[,common_cell]
  cell.embeddings_sub = cell.embeddings[common_cell,]
  cell.embeddings_sub = rownames_to_column(as.data.frame(cell.embeddings_sub), var = 'cell_id')
  
  exp = rownames_to_column(exp, var = 'tr_id')
  
  exp_df = exp %>% pivot_longer(cols = -tr_id, names_to = 'cell_id', values_to = 'exp') %>%
    pivot_wider(names_from = tr_id, values_from = exp)
  
  cell.embeddings_sub = left_join(cell.embeddings_sub, exp_df, by = "cell_id")
  cell.embeddings_sub = na.omit(cell.embeddings_sub)

  cell.embeddings_sub[[ge]] = rowSums(cell.embeddings_sub[, trs, drop = FALSE], na.rm = TRUE)
  p1 = ggplot( cell.embeddings_sub %>%
                 arrange(ge), aes(UMAP_1, UMAP_2, color = get(ge))) +
    geom_point() +
    scale_color_gradient(low = "gray", high = "purple") +
    theme_minimal() +
    labs(color = paste0(ge, ' Expression'))
  
  ##
  
  selected_cells = colnames(exp)[-1][apply(exp[, -1], 2, function(x) sum(x > 0) > 0)]
  
  cell.embeddings_sub = as.data.frame(as.data.frame(cell.embeddings)[selected_cells, ])
  cell.embeddings_sub$cell_id = rownames(cell.embeddings_sub)
  exp = exp[, c("tr_id", selected_cells)]
  exp_df = exp %>% pivot_longer(cols = -tr_id, names_to = 'cell_id', values_to = 'exp') %>%
    pivot_wider(names_from = tr_id, values_from = exp)
  
  cell.embeddings_sub = left_join(cell.embeddings_sub, exp_df, by = "cell_id")
  

  cell.embeddings_sub$class = apply(cell.embeddings_sub[, trs, drop = FALSE], 1, function(x) {
    if (sum(x > 0) == 1) {
      paste0('only_', names(which(x > 0)))
    } else if (sum(x > 0) > 1) {
      'multiple_exp'
    } else {
      'no_exp'
    }
  })
  
  p2 = ggplot(cell.embeddings_sub, aes(UMAP_1, UMAP_2, color = class)) +
    geom_point(alpha = 0.7) +
    theme_minimal()
  
  class_counts <- cell.embeddings_sub %>%
    count(class) %>%
    mutate(percentage = n / sum(n) * 100,
           label = paste0(class, " (", n, " - ", round(percentage, 1), "%)"))

  p3 <- ggplot(class_counts, aes(x = "", y = n, fill = label)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y") +
    theme_void() +
    theme(legend.position = "right") +
    labs(fill = "Class")

  
  print(p / p1 / p2 / p3)
}
dev.off()


three_tr_ges = names(table(id_map$gene_id)[table(id_map$gene_id) == 3])
three_tr_ges = id_map[id_map$gene_id %in% three_tr_ges,]
cells_stat_df <- data.frame()


for (ge in unique(three_tr_ges$gene_id)) {
  trs <- three_tr_ges[three_tr_ges$gene_id == ge,]$transcript_id
  exp <- as.data.frame(sce.nano@assays$RNA@counts[trs, ])

  cells_stat <- lapply(trs, function(tr) {
    colnames(exp)[exp[tr, ] > 0]
  })
  names(cells_stat) <- trs

  unique_cells_stat <- lapply(trs, function(tr) {
    setdiff(cells_stat[[tr]], unlist(cells_stat[setdiff(trs, tr)]))
  })
  names(unique_cells_stat) <- trs
  

  all_cells <- unique(unlist(cells_stat))
  all_cell_num <- length(all_cells)
  
  cells_stat_sub <- data.frame(
    gene_id = ge,
    all_cell_num = all_cell_num,
    tr_1_id = trs[1],
    tr_1_num = length(unique_cells_stat[[1]]),
    tr_1_ratio = length(unique_cells_stat[[1]]) / all_cell_num,
    tr_2_id = trs[2],
    tr_2_num = length(unique_cells_stat[[2]]),
    tr_2_ratio = length(unique_cells_stat[[2]]) / all_cell_num,
    tr_3_id = trs[3],
    tr_3_num = length(unique_cells_stat[[3]]),
    tr_3_ratio = length(unique_cells_stat[[3]]) / all_cell_num
  )
  
  cells_stat_df <- rbind(cells_stat_df, cells_stat_sub)
}

write.table(cells_stat_df,'three_tr_genes.txt',sep = '\t',row.names = F,quote = F)
#---------------------------------------------------------------------------------------------------------------

  



# three example
cells_stat_df = read.table('three_tr_genes.txt',sep = '\t',header = T)
cells_stat_df_select = cells_stat_df[cells_stat_df$tr_1_ratio > 0.2& cells_stat_df$tr_2_ratio > 0.2 & cells_stat_df$tr_3_ratio > 0.2,]
cells_stat_df_select = left_join(cells_stat_df_select,unique(id_map[,c('gene_id','gene_type')]))
genes = cells_stat_df_select[cells_stat_df_select$gene_type == 'Protein coding',]$gene_id
genes = c('ENSG00000178074.6','ENSG00000122877.17','ENSG00000196110.8')

cell.embeddings = sce.ill@reductions[["umap"]]@cell.embeddings

pdf('three_example_featureplot.pdf',width = 15,height = 10)
for (ge in genes) {
  trs = id_map[id_map$gene_id %in% ge, ]$transcript_id
  ge_exons = dplyr::filter(gtf, transcript_id %in% trs &type == "exon")
  ge_exons$transcript_id = factor(ge_exons$transcript_id,levels = trs)
  p = ggplot(ge_exons,aes(xstart = start,xend = end,y = transcript_id
  )) +
    geom_range()+
    geom_intron(
      data = to_intron(ge_exons, "transcript_id"),
      aes(strand = strand)
    )+
    theme_bw()
  exp = as.data.frame(sce.nano@assays$RNA@data[trs, ])
  common_cell = intersect(rownames(cell.embeddings),colnames(exp))
  exp = exp[,common_cell]
  cell.embeddings_sub = cell.embeddings[common_cell,]
  cell.embeddings_sub = rownames_to_column(as.data.frame(cell.embeddings_sub), var = 'cell_id')
  
  exp = rownames_to_column(exp, var = 'tr_id')
  
  exp_df = exp %>% pivot_longer(cols = -tr_id, names_to = 'cell_id', values_to = 'exp') %>%
    pivot_wider(names_from = tr_id, values_from = exp)
  
  cell.embeddings_sub = left_join(cell.embeddings_sub, exp_df, by = "cell_id")
  cell.embeddings_sub = na.omit(cell.embeddings_sub)

  cell.embeddings_sub[[ge]] = rowSums(cell.embeddings_sub[, trs, drop = FALSE], na.rm = TRUE)
  p1 = ggplot( cell.embeddings_sub %>%
                 arrange(ge), aes(UMAP_1, UMAP_2, color = get(ge))) +
    geom_point() +
    scale_color_gradient(low = "gray", high = "purple") +
    theme_minimal() +
    labs(color = paste0(ge, ' Expression'))
  
  ##
  
  selected_cells = colnames(exp)[-1][apply(exp[, -1], 2, function(x) sum(x > 0) > 0)]
  
  cell.embeddings_sub = as.data.frame(as.data.frame(cell.embeddings)[selected_cells, ])
  cell.embeddings_sub$cell_id = rownames(cell.embeddings_sub)
  exp = exp[, c("tr_id", selected_cells)]
  exp_df = exp %>% pivot_longer(cols = -tr_id, names_to = 'cell_id', values_to = 'exp') %>%
    pivot_wider(names_from = tr_id, values_from = exp)
  
  cell.embeddings_sub = left_join(cell.embeddings_sub, exp_df, by = "cell_id")
  

  cell.embeddings_sub$class = apply(cell.embeddings_sub[, trs, drop = FALSE], 1, function(x) {
    if (sum(x > 0) == 1) {
      paste0('only_', names(which(x > 0)))
    } else if (sum(x > 0) > 1) {
      'multiple_exp'
    } else {
      'no_exp'
    }
  })
  
  p2 = ggplot(cell.embeddings_sub, aes(UMAP_1, UMAP_2, color = class)) +
    geom_point(alpha = 0.7) +
    theme_minimal()
  
  class_counts <- cell.embeddings_sub %>%
    count(class) %>%
    mutate(percentage = n / sum(n) * 100,
           label = paste0(class, " (", n, " - ", round(percentage, 1), "%)"))
  

  p3 <- ggplot(class_counts, aes(x = "", y = n, fill = label)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y") +
    theme_void() +
    theme(legend.position = "right") +
    labs(fill = "Class")
  
  
  print(p / p1 / p2 / p3)
}
dev.off()
