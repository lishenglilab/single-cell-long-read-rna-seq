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
