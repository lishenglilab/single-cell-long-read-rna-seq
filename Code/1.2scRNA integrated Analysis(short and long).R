library(tidyverse)
library(Seurat)
library(harmony)

#short read
dir1 <- "Data/short/single_cell/"
## get sample list
samples = list.files(dir1)

#create Seurat objects
samList = lapply(samples, function(sp){
  folder=file.path(dir1, sp)
  CreateSeuratObject(counts = Read10X(folder),
                     project = sp,
                     min.cells = 10,
                     min.features = 200)
})
##filter out doublets
# run 1.Scrublet.py
dir_scrublet="Result/Scrublet/" # set estimated doublet rate at 6 %
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

## filter cells
rt_summary <- data.frame()
names(samList)=samples
for(i in samples){
  samList[[i]] <- PercentageFeatureSet(object = samList[[i]],  pattern = "^MT-", col.name = "percent.mt") #rat mitochondrial genes start with Mt- 
  cell_before <- nrow(samList[[i]]@meta.data)
  scrublet_result=read.table(file=paste0(dir_scrublet, "/", i,"doublets.txt"),header=F,sep="\t",check.names=F)
  colnames(scrublet_result)=c("barcode","double cell")
  scrublet_result=scrublet_result[match(rownames(samList[[i]]@meta.data), scrublet_result$barcode),]
  samList[[i]]@meta.data$double.cell=scrublet_result$`double cell`
  
  # doublet rate
  if("True" %in% unique(samList[[i]]@meta.data$double.cell)){
    doublet_rate <- round(table(samList[[i]]@meta.data$double.cell)["True"]*100/cell_before, 2)
  }else{
    doublet_rate = 0
  }
  
  nCount=samList[[i]]@meta.data$nCount_RNA
  nFeature=samList[[i]]@meta.data$nFeature_RNA
  mt=samList[[i]]@meta.data$percent.mt
  samList[[i]] <- AddMetaData(samList[[i]], log10(nFeature), col.name = "log.nFeature")
  samList[[i]] <- AddMetaData(samList[[i]], log10(nCount), col.name = "log.nCount_RNA")
  samList[[i]] <- subset(samList[[i]],
                         subset =
                           log.nFeature > median(log10(nFeature))-3*mad(log10(nFeature)) &
                           log.nCount_RNA > median(log10(nCount))-3*mad(log10(nCount)) &
                           percent.mt < median(mt) + 3*mad(mt) &
                           double.cell !="True")
  cell_after <- nrow(samList[[i]]@meta.data)
  
  # summary
  rt_summary_i <- data.frame(Sample = i, Before_Filter = cell_before, After_Filter = cell_after, Doublet_Rate = doublet_rate, Cell_Reduction_Rate = (cell_before-cell_after)*100/cell_before)
  rt_summary <- rbind(rt_summary, rt_summary_i)
  
  cat(i, "Before filter :", cell_before, "cells; ", "After filter :", cell_after,"cells\n",
      "Doublet rate: ", doublet_rate, "%; ", "Cell reduction rate: ", (cell_before-cell_after)*100/cell_before, "%\n")
  
}

samList <- lapply(X = samList, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = F) # Assign Cell-Cycle Scores
})
features <- SelectIntegrationFeatures(object.list = samList) #set lower features to avoid integration failure
samList <- lapply(X = samList, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE) 
  x <- RunPCA(x, features = features, verbose = FALSE)
})


#### integrate samples ####
samList.anchors <- FindIntegrationAnchors(object.list = samList, anchor.features = features, reduction = "rpca")
sce.mer <- IntegrateData(anchorset = samList.anchors)

## cluster
DefaultAssay(sce.mer) <- "integrated"
sce.mer <- ScaleData(sce.mer, verbose=FALSE) # normally do not need regress out the sources of heterogeneity from the data
sce.mer <- RunPCA(sce.mer, npcs=30,verbose=FALSE)
sce.mer <- RunUMAP(sce.mer,reduction="pca",dims=1:30)
sce.mer <- FindNeighbors(sce.mer, reduction="pca",dims = 1:30)
sce.mer <- FindClusters(sce.mer, resolution = 0.8)
write.table(rt_summary, file = "Result/filtering_summary.txt", sep = "\t", row.names = F, quote = F)
saveRDS(sce.mer, file="Result/integration_Illumina.rds")


## Automatic annotation
library(SingleR)
refdata <- get(load("Reference/ref_Human_all.RData"))
testdata <- GetAssayData(sce.mer, slot="data")
clusters <- sce.mer@meta.data$seurat_clusters
cellpred <- SingleR(test = testdata, ref = refdata, 
                    labels = refdata$label.fine, clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")

write.table(celltype,file="Result/singleR-anotation_sum.txt",sep='\t',row.names=F,quote = F)


## long read 
data_dir <- "Data/long/single_cell/"
samples <- c("AD1","N1","AD2","N2","AD3","N3")
samList = lapply(samples, function(sp){
  folder=paste0(data_dir, sp,'/FLAMES_out/transcript_count.csv.gz')
  count <- data.table::fread(folder,sep = ",",check.names = F,data.table = F)
  count$gene_id <- NULL
  rownames(count) <- count$transcript_id
  count$transcript_id <- NULL
  colnames(count) <- paste0(colnames(count),'-1')
  CreateSeuratObject(counts = count,
                     project = sp)
})

samList <- lapply(X = samList, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = samList) #set lower features to avoid integration failure
samList <- lapply(X = samList, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE) 
  x <- RunPCA(x, features = features, verbose = FALSE)
})


#### integrate samples ####
samList.anchors <- FindIntegrationAnchors(object.list = samList, anchor.features = features, reduction = "rpca")
sce.mer <- IntegrateData(anchorset = samList.anchors)
## cluster
DefaultAssay(sce.mer) <- "integrated"
sce.mer <- ScaleData(sce.mer, verbose=FALSE) # normally do not need regress out the sources of heterogeneity from the data
sce.mer <- RunPCA(sce.mer, npcs=30,verbose=FALSE)
sce.mer <- RunUMAP(sce.mer,reduction="pca",dims=1:30)
sce.mer <- FindNeighbors(sce.mer, reduction="pca",dims = 1:30)
sce.mer <- FindClusters(sce.mer, resolution = 0.8)
saveRDS(sce.mer, file="Result/integration_nanopore.rds")


# correlations
nano_exp = sce_nano@assays$RNA@data
ill_exp = sce_ill@assays$RNA@data
common_barcodes = intersect(colnames(nano_exp), colnames(ill_exp))
nano_exp = nano_exp[, common_barcodes]
ill_exp = ill_exp[, common_barcodes]
nano_exp$transcript_id = rownames(nano_exp)
nano_exp <- nano_exp %>%
mutate(gene_id = ifelse(str_starts(transcript_id, 'ENSG'),
gsub('-.*', '', transcript_id),
id_map$gene_id[match(transcript_id, id_map$transcript_id)]))
nano_exp = left_join(nano_exp,unique(id_map[,c('gene_id','gene_name')]))
nano_exp$gene_id = NULL
nano_exp$transcript_id = NULL
common_genes = intersect(rownames(nano_exp_sum), rownames(ill_exp))
nano_exp = nano_exp[nano_exp$gene_name %in% common_genes,]

nano_exp_sum = fread('nano_ge_exp_summed.csv',sep = ',',data.table = F)
rownames(nano_exp_sum) = nano_exp_sum$gene_name
nano_exp_sum$gene_name = NULL
ill_exp = ill_exp[common_genes,]
nano_exp_sum = nano_exp_sum[common_genes,]
ill_exp = as.data.frame(as.matrix(ill_exp))
correlations <- sapply(common_barcodes, function(barcode) {
  cor(nano_exp_sum[, barcode], ill_exp[, barcode], method = "pearson")
})

median_cor <- median(correlations)
set.seed(123)  
random_correlations <- sapply(1:length(common_barcodes), function(i) {
  random_barcode <- sample(colnames(ill_exp), 1)
  cor(nano_exp_sum[, common_barcodes[i]], ill_exp[, random_barcode], method = "pearson")
})


median_random_cor <- median(random_correlations)

plot_data <- data.frame(
  Type = rep(c("Matched", "Random"), each = length(correlations)),
  Correlation = c(correlations, random_correlations)
)
correlation <- cor(gene_per_cell_sub$Nano, gene_per_cell_sub$Ill)
p = ggplot(plot_data, aes(x = Type, y = Correlation, fill = Type)) +
  stat_boxplot(geom = 'errorbar',width = 0.2)+
  geom_boxplot(outlier.colour = NA) +
  stat_summary(fun = median, geom = "point", shape = 23, size = 4, fill = "red") +
  theme_minimal() +
  labs(x = "Cell Pairing Type",
       y = "Pearson Correlation")+
  stat_compare_means(comparisons = list(c('Matched','Random')))
ggsave(plot = p,'ALL_exp_pearson_cor_boxplot.pdf',width = 6,height = 4)

merged_gene_per_cell = read.table('/data/haowu/sc_long/nano_ill_cor/merged_gene_per_cell.txt',sep = '\t',header = T)
correlation <- cor(merged_gene_per_cell$Nano, merged_gene_per_cell$Ill)
p = ggplot(merged_gene_per_cell, aes(x = Nano, y = Ill)) +
  geom_point(alpha = 0.8) +
  geom_text(x = min(merged_gene_per_cell$Nano), y = max(merged_gene_per_cell$Ill), label = paste("Pearson correlation =", round(correlation, 2)), hjust = 0, vjust = 1)+
  theme_bw()
ggsave(plot = p,'ALL_gene_per_cell_cor_scatter.pdf',width = 6,height = 4)
