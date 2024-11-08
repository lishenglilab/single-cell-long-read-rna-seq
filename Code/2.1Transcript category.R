rm(list = ls()) # remove environment variables
library(tidyverse)
library(rtracklayer)

# set working direction
setwd("Result/CAGE_chromatin/")
## extract single cell transcript TSS range bed
# import single cell object
sce.mer <- read_rds("Result/integration_nanopore.rds")

# import gtf data
gtf <- import("Result/isoform_annotated.merged.sorted.gtf")
gtf$transcript_id <- gsub("_","-",gtf$transcript_id)
gtf <- gtf[gtf$transcript_id %in% rownames(sce.mer@assays$RNA@data)]
length(unique(gtf$transcript_id))

# extract transcript information
index <- which(gtf$type == "transcript") # extract transcript index
tr_info <- data.frame(Gene_id = gtf$gene_id[index],
                   Transcript_Id = gtf$transcript_id[index],
                   gtf@seqnames[index],
                   gtf@ranges[index],
                   gtf@strand[index])
colnames(tr_info)[3:7] <- c("Chr","Start","End","Width","Strand")

# extract transcript TSS
transcript_positive_TSS <- tr_info[tr_info$Strand == "+",]
transcript_positive_TSS_bed <- data.frame(transcript_positive_TSS$Chr, 
                                          transcript_positive_TSS$Start - 500,
                                          transcript_positive_TSS$Start + 500,
                                          transcript_positive_TSS$Transcript_Id,
                                          0,
                                          transcript_positive_TSS$Strand)
transcript_positive_TSS_bed$transcript_positive_TSS.Start...500[transcript_positive_TSS_bed$transcript_positive_TSS.Start...500 < 0] <- 1

transcript_negative_TSS <- tr_info[tr_info$Strand == "-",]
transcript_negative_TSS$TSS <- transcript_negative_TSS$End
transcript_negative_TSS_bed <- data.frame(transcript_negative_TSS$Chr, 
                                          transcript_negative_TSS$End - 500,
                                          transcript_negative_TSS$End + 500,
                                          transcript_negative_TSS$Transcript_Id,
                                          0,
                                          transcript_negative_TSS$Strand)

# merge data
colnames(transcript_positive_TSS_bed) <- ""
colnames(transcript_negative_TSS_bed) <- ""
transcript_TSS_bed <- rbind(transcript_positive_TSS_bed, transcript_negative_TSS_bed)

# export bed
write.table(transcript_TSS_bed, "isoform_annotated.merged_transcript_TSS.bed", sep = "\t", col.names = F, row.names = F, quote = F)

## extract genecode v38 transcript TSS range bed
# import gtf data
gtf <- import("~/reference/Gencode/Human/Gencode_v38/gencode.v38.primary_assembly.annotation.gtf")
gtf <- gtf[gtf$transcript_id %in% rownames(sce.mer@assays$RNA@data)]
length(unique(gtf$transcript_id))

# extract transcript information
index <- which(gtf$type == "transcript") # extract transcript index
tr_info <- data.frame(Gene_id = gtf$gene_id[index],
                      Transcript_Id = gtf$transcript_id[index],
                      gtf@seqnames[index],
                      gtf@ranges[index],
                      gtf@strand[index])
colnames(tr_info)[3:7] <- c("Chr","Start","End","Width","Strand")

# extract transcript TSS
transcript_positive_TSS <- tr_info[tr_info$Strand == "+",]
transcript_positive_TSS_bed <- data.frame(transcript_positive_TSS$Chr, 
                                          transcript_positive_TSS$Start - 500,
                                          transcript_positive_TSS$Start + 500,
                                          transcript_positive_TSS$Transcript_Id,
                                          0,
                                          transcript_positive_TSS$Strand)
transcript_positive_TSS_bed$transcript_positive_TSS.Start...500[transcript_positive_TSS_bed$transcript_positive_TSS.Start...500 < 0] <- 1

transcript_negative_TSS <- tr_info[tr_info$Strand == "-",]
transcript_negative_TSS$TSS <- transcript_negative_TSS$End
transcript_negative_TSS_bed <- data.frame(transcript_negative_TSS$Chr, 
                                          transcript_negative_TSS$End - 500,
                                          transcript_negative_TSS$End + 500,
                                          transcript_negative_TSS$Transcript_Id,
                                          0,
                                          transcript_negative_TSS$Strand)


# merge data
colnames(transcript_positive_TSS_bed) <- ""
colnames(transcript_negative_TSS_bed) <- ""
transcript_TSS_bed <- rbind(transcript_positive_TSS_bed, transcript_negative_TSS_bed)

# export bed
write.table(transcript_TSS_bed, "gencode_v38_transcript_TSS.bed", sep = "\t", col.names = F, row.names = F, quote = F)

## extract single cell gene TSS range bed
# import gtf data
gtf <- import("Result/isoform_annotated.merged.sorted.gtf")
gtf$transcript_id <- gsub("_","-",gtf$transcript_id)
gtf <- gtf[gtf$transcript_id %in% rownames(sce.mer@assays$RNA@data)]
length(unique(gtf$transcript_id))

# extract transcript information
index <- which(gtf$type == "transcript") # extract transcript index
tr_info <- data.frame(Gene_id = gtf$gene_id[index],
                      Transcript_Id = gtf$transcript_id[index],
                      gtf@seqnames[index],
                      gtf@ranges[index],
                      gtf@strand[index])
colnames(tr_info)[3:7] <- c("Chr","Start","End","Width","Strand")

# extract gene TSS
gene_positive_TSS <- tr_info[tr_info$Strand == "+",] %>% 
  group_by(Gene_id) %>% 
  arrange(Start) %>% 
  summarise(TSS = .data$Start[1])

gene_negative_TSS <- tr_info[tr_info$Strand == "-",] %>% 
  group_by(Gene_id) %>% 
  arrange(desc(End)) %>% 
  summarise(TSS = .data$End[1])

double_strand_gene = intersect(gene_positive_TSS$Gene_id, gene_negative_TSS$Gene_id) # check gene strand
gene_positive_TSS <- gene_positive_TSS[!(gene_positive_TSS$Gene_id %in% double_strand_gene),]
gene_negative_TSS <- gene_negative_TSS[!(gene_negative_TSS$Gene_id %in% double_strand_gene),]

# build gene bed
gene_positive_TSS <- left_join(gene_positive_TSS, unique(tr_info[, c(1, 3, 7)]), by = "Gene_id")
gene_positive_TSS <- gene_positive_TSS[gene_positive_TSS$Strand != "-",]
gene_positive_TSS_bed <- data.frame(gene_positive_TSS$Chr, 
                                    gene_positive_TSS$TSS - 500, 
                                    gene_positive_TSS$TSS + 500, 
                                    gene_positive_TSS$Gene_id, 
                                    0, 
                                    gene_positive_TSS$Strand)
gene_positive_TSS_bed$gene_positive_TSS.TSS...500[gene_positive_TSS_bed$gene_positive_TSS.TSS...500 < 0] <- 1

gene_negative_TSS <- left_join(gene_negative_TSS, unique(tr_info[, c(1, 3, 7)]), by = "Gene_id")
gene_negative_TSS <- gene_negative_TSS[gene_negative_TSS$Strand != "+",]
gene_negative_TSS_bed <- data.frame(gene_negative_TSS$Chr, 
                                    gene_negative_TSS$TSS - 500, 
                                    gene_negative_TSS$TSS + 500, 
                                    gene_negative_TSS$Gene_id, 
                                    0, 
                                    gene_negative_TSS$Strand)
# merge data
colnames(gene_positive_TSS_bed) <- ""
colnames(gene_negative_TSS_bed) <- ""
gene_TSS_bed <- rbind(gene_positive_TSS_bed, gene_negative_TSS_bed)

# export bed
write.table(gene_TSS_bed, "sc_gene_TSS.bed", sep = "\t", col.names = F, row.names = F, quote = F)

##
cd Result/CAGE_chromatin
# Roadmap intersection
#bedtools intersect -a sc_gene_TSS.bed -b Roadmap_chromatin_state.bed -wo > sc_gene_TSS_Roadmap_intersect

# CAGE intersection
#bedtools intersect -a gencode_v38_transcript_TSS.bed -b hg38.cage_peak.bed -s -wo > gencode_v38_transcript_TSS_CAGE_intersect
#bedtools intersect -a isoform_annotated.merged_transcript_TSS.bed -b hg38.cage_peak.bed -s -wo > isoform_annotated.merged_transcript_TSS_CAGE_intersect

rm(list = ls()) # remove environment variables
library(tidyverse)

# import transcript intersect data
sc_transcript_cage <- data.table::fread("isoform_annotated.merged_transcript_TSS_CAGE_intersect", header = F, data.table = F, check.names = F)
gencode_transcript_cage <- data.table::fread("gencode_v38_transcript_TSS_CAGE_intersect", header = F, data.table = F, check.names = F)

# import gene intersect data
gene_roadmap <- data.table::fread("sc_gene_TSS_Roadmap_intersect", header = F, data.table = F, check.names = F)

## CAGE intersection
transcript_cage <- rbind(sc_transcript_cage, gencode_transcript_cage)
transcript_cage$id <- paste0(transcript_cage$V7, ":", transcript_cage$V8, "|", transcript_cage$V9, ":", transcript_cage$V12)
transcript_cage <- unique(transcript_cage)
names(which(table(paste(transcript_cage$V4, transcript_cage$id)) > 1))

transcript_cage_summary <- transcript_cage %>% 
  group_by(V4) %>% 
  summarise(Intersected_CAGE_peaks = paste(id, collapse = "; "))
colnames(transcript_cage_summary)[1] <- "Transcript_id"

# export data
write.table(transcript_cage_summary, "transcript_CAGE_intersect_summary.txt", sep = "\t", row.names = F, quote = F)

## gene chromatin state intersection
gene_roadmap_summary <- gene_roadmap %>% 
  group_by(V4, V6) %>% 
  summarise(Intersected_chromatin_state = paste(unique(V10), collapse = "; "))
colnames(gene_roadmap_summary)[1:2] <- c("Gene_id", "Strand")
names(which(table(gene_roadmap_summary$Gene_id) > 1))

transcript_info <- data.table::fread("~/project/single_cell_full_length/Result/id_map.txt", header = T, data.table = F, check.names = F)
setdiff(gene_roadmap_summary$Gene_id, transcript_info$gene_id)
gene_roadmap_summary <- left_join(gene_roadmap_summary, transcript_info[, 1:2], by = c("Gene_id" = "gene_id")) %>% 
  select(Gene_id, transcript_id, Strand, Intersected_chromatin_state)
table(is.na(gene_roadmap_summary$transcript_id))

# export data
write.table(gene_roadmap_summary, "gene_chromatin_state_intersect_summary.txt", sep = "\t", row.names = F, quote = F)

transcript_cage_summary <- data.table::fread("transcript_CAGE_intersect_summary.txt", header = T, data.table = F, check.names = F)
names(which(table(transcript_cage_summary$Transcript_id) > 1))

gene_roadmap_summary <- data.table::fread("gene_chromatin_state_intersect_summary.txt", header = T, data.table = F, check.names = F)
names(which(table(gene_roadmap_summary$transcript_id) > 1))

transcript_info <- data.table::fread("Result/id_map.txt", header = T, data.table = F, check.names = F)
transcript_info$Annotation_Status <- ifelse(str_detect(transcript_info$transcript_id,pattern = "-"),"unannotated","annotated")
## transcript CAGE + gene chromatin state
RNA_cage_DNA_roadmap <- data.frame(Transcript_id = transcript_info$transcript_id[transcript_info$Annotation_Status == "unannotated"])
RNA_cage_DNA_roadmap <- data.frame(Transcript_id = RNA_cage_DNA_roadmap[RNA_cage_DNA_roadmap$Transcript_id %in% rownames(sce.nano@assays$RNA@data),])
RNA_cage_DNA_roadmap <- left_join(RNA_cage_DNA_roadmap, transcript_cage_summary, by = "Transcript_id")
RNA_cage_DNA_roadmap <- left_join(RNA_cage_DNA_roadmap, gene_roadmap_summary[c(2, 4)], by = c("Transcript_id" = "transcript_id"))
names(which(table(RNA_cage_DNA_roadmap$Transcript_id) > 1))

RNA_cage_DNA_roadmap$CAGE_status <- ifelse(!is.na(RNA_cage_DNA_roadmap$Intersected_CAGE_peaks), "CAGE", "")
RNA_cage_DNA_roadmap$chromatin_status <- ifelse(!is.na(RNA_cage_DNA_roadmap$Intersected_chromatin_state), "Transcript_activation", "")
RNA_cage_DNA_roadmap$Status <- paste(RNA_cage_DNA_roadmap$CAGE_status, RNA_cage_DNA_roadmap$chromatin_status)
RNA_cage_DNA_roadmap$Status[RNA_cage_DNA_roadmap$Status == " "] <- "None"
table(RNA_cage_DNA_roadmap$Status)
RNA_cage_DNA_roadmap <- RNA_cage_DNA_roadmap[!duplicated(RNA_cage_DNA_roadmap[,c("Transcript_id","Status")]),]
# plot pie chart
summary_pie <- data.frame(table(RNA_cage_DNA_roadmap$Status))
colnames(summary_pie)[1] <- "transcription_status"
summary_pie$percentage <- summary_pie$Freq*100/sum(summary_pie$Freq)
myPalette <- brewer.pal(11, "Set3")
pdf("novel_transcription_status_DNA_pie.pdf", width = 12, height = 8)
pie(summary_pie$Freq , labels = paste0(summary_pie$transcription_status, "(", summary_pie$Freq, " ", round(summary_pie$percentage, 2), "%)"), border = "white", col = myPalette )
dev.off()

## annotated transcript CAGE + gene chromatin state
annotated_RNA_cage_DNA_roadmap <- data.frame(Transcript_id = transcript_info$transcript_id[transcript_info$Annotation_Status == "annotated"])
annotated_RNA_cage_DNA_roadmap <- data.frame(Transcript_id = annotated_RNA_cage_DNA_roadmap[annotated_RNA_cage_DNA_roadmap$Transcript_id %in% rownames(sce.nano@assays$RNA@data),])
annotated_RNA_cage_DNA_roadmap <- left_join(annotated_RNA_cage_DNA_roadmap, transcript_cage_summary, by = "Transcript_id")
annotated_RNA_cage_DNA_roadmap <- left_join(annotated_RNA_cage_DNA_roadmap, gene_roadmap_summary[c(2, 4)], by = c("Transcript_id" = "transcript_id"))
names(which(table(annotated_RNA_cage_DNA_roadmap$Transcript_id) > 1))

annotated_RNA_cage_DNA_roadmap$CAGE_status <- ifelse(!is.na(annotated_RNA_cage_DNA_roadmap$Intersected_CAGE_peaks), "CAGE", "")
annotated_RNA_cage_DNA_roadmap$chromatin_status <- ifelse(!is.na(annotated_RNA_cage_DNA_roadmap$Intersected_chromatin_state), "Transcript_activation", "")
annotated_RNA_cage_DNA_roadmap$Status <- paste(annotated_RNA_cage_DNA_roadmap$CAGE_status, annotated_RNA_cage_DNA_roadmap$chromatin_status)
annotated_RNA_cage_DNA_roadmap$Status[annotated_RNA_cage_DNA_roadmap$Status == " "] <- "None"
table(annotated_RNA_cage_DNA_roadmap$Status)
annotated_RNA_cage_DNA_roadmap <- annotated_RNA_cage_DNA_roadmap[!duplicated(annotated_RNA_cage_DNA_roadmap[,c("Transcript_id","Status")]),]
# plot pie chart
annotated_summary_pie <- data.frame(table(annotated_RNA_cage_DNA_roadmap$Status))
colnames(annotated_summary_pie)[1] <- "transcription_status"
annotated_summary_pie$percentage <- annotated_summary_pie$Freq*100/sum(annotated_summary_pie$Freq)

myPalette <- brewer.pal(11, "Set3")
pdf("annotated_transcription_status_DNA_pie.pdf", width = 12, height = 8)
pie(annotated_summary_pie$Freq , labels = paste0(annotated_summary_pie$transcription_status, "(", annotated_summary_pie$Freq, " ", round(annotated_summary_pie$percentage, 2), "%)"), border = "white", col = myPalette )
dev.off()

# fisher test
x <- cbind(summary_pie, annotated_summary_pie)
rownames(x) <- x[, 1]
x <- x[, c(2, 5)]
colnames(x) <- c("unannotated", "annotated")

chisq <- chisq.test(x)
print(x)
print(paste(chisq$statistic, chisq$p.value))


#######
transcript_assembly_classification <- read.table("Result/new_AD_classification.txt",sep = "\t",header = T)
transcript_assembly_classification <- transcript_assembly_classification[transcript_assembly_classification$isoform %in% rownames(sce.nano),]
transcript_assembly_classification$structural_category [!(transcript_assembly_classification$structural_category %in% c("full-splice_match","incomplete-splice_match","novel_in_catalog","novel_not_in_catalog"))] <- "Other"
res <- data.frame()
for (type in unique(transcript_assembly_classification$structural_category)) {
  df <- data.frame(Transcript_id = transcript_assembly_classification[transcript_assembly_classification$structural_category == type,]$isoform)
  df <- data.frame(Transcript_id = df[df$Transcript_id %in% rownames(sce.nano@assays$RNA@data),])
  df <- left_join(df, transcript_cage_summary, by = "Transcript_id")
  df <- left_join(df, gene_roadmap_summary[c(2, 4)], by = c("Transcript_id" = "transcript_id"))
  names(which(table(df$Transcript_id) > 1))
  
  df$CAGE_status <- ifelse(!is.na(df$Intersected_CAGE_peaks), "CAGE", "")
  df$chromatin_status <- ifelse(!is.na(df$Intersected_chromatin_state), "Transcript_activation", "")
  df$Status <- paste(df$CAGE_status, df$chromatin_status)
  df$Status[df$Status == " "] <- "None"
  table(df$Status)
  df <- df[!duplicated(df[,c("Transcript_id","Status")]),]
  res_sub <- df %>% group_by(Status) %>% summarise(Num = n())
  res_sub$type <- type
  res <- rbind(res,res_sub)
}
ggplot(res,aes(type,Num,fill = Status))+
  geom_bar(stat = "identity")+
  theme_classic()+
  coord_flip()
ggsave("categories_transcription_status_DNA_bar_count.pdf",width = 10,height = 4)
p =ggplot(res,aes(type,Num,fill = Status))+
  geom_bar(stat = "identity",position = "fill")+
  theme_classic()+
  coord_flip()
ggsave("categories_transcription_status_DNA_bar_pct.pdf",width = 10,height = 4)

data = p$data
total_by_type <- data %>%
  group_by(type) %>%
  summarise(total = sum(Num))

percentage_data <- data %>%
  group_by(type, Status) %>%
  summarise(count = sum(Num)) %>%
  left_join(total_by_type, by = "type") %>%
  mutate(percentage = count / total * 100)
write.table(percentage_data,'Reslt/categories_transcription_status_polt_data.txt',sep = '\t',row.names = F,quote = F)

