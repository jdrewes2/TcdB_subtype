---
title: "DESEQ2 analysis"
author: "Julia Drewes"
date: "2024-11-07"
operating system: "MacBook Pro Sonoma 14.6.1"
R version: "4.4.1"
output: html_document
---
  
  
#Colon data normalized across all samples
  ```{r}
rm(list=ls())
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("DESeq2")
library(readxl)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(dplyr)

group = as.data.frame(read_excel("./metadata.xlsx", sheet = 1))

rownames(group)= paste0(group$`Novogene sampleID`,"")

rna_count = read.csv("./gene_count.csv",check.names = F,row.names = 1)
# extract information of geneid and genename
info = as.data.frame(cbind(gene_id=rownames(rna_count),gene_name=rna_count$gene_name))
rownames(info)=info$gene_id

head(rna_count)

rna_input=rna_count[,rownames(group)]

#colData = group --> Data frame group
#design = ~ group --> column group

dds <- DESeqDataSetFromMatrix(countData = rna_input,
                              colData = group,
                              design = ~ group)
dds
# fitering > 10 and smallest samplesize per group (here we use 6)
keep <- rowSums(counts(dds) >= 10) >=6

dds <- dds[keep,]

# Differential expression analysis
dds <- DESeq(dds)
# to retrieve the normalized counts matrix from dds, we use the counts() function and add the argument normalized=TRUE
normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts=as.data.frame(normalized_counts)
## map gene name to gene id
info_tmp=info[rownames(normalized_counts),]
normalized_counts$gene_name = info_tmp$gene_name
##exporting output DESEQ2 table to Desktop
write.csv(normalized_counts, "./Normalized_raw_gene-counts.csv")

######## Pairwise comparisons
# 1 comparison: TcdB1vsSham = get_pairwise
result1 = results(dds, contrast = c("group","TcdB1","Sham"))
output1= as.data.frame(result1)
# map gene name to gene id
info_tmp=info[rownames(output1),]
output1$gene_name = info_tmp$gene_name
#exporting log2Foldchange table
write.csv(output1, "./log2Foldchange_padj-table/TcdB1_vs_Sham.csv")

# sig gene: abs(log2fc) > 1 && padj<0.05
output1$diffexpressed <- NA
output1$diffexpressed[output1$log2FoldChange > 1 & output1$padj < 0.05] <- "UP"
output1$diffexpressed[output1$log2FoldChange < -1 & output1$padj < 0.05] <- "DOWN"
mycolors <- c("blue", "red", "grey")
names(mycolors) <- c("DOWN", "UP", NA)
output1$delabel <- NA
output1$delabel[!is.na(output1$diffexpressed)] <- output1$gene_name[!is.na(output1$diffexpressed)]
p <- ggplot(data=output1, aes(x=log2FoldChange, y=-log10(padj),color=diffexpressed,label=delabel)) + geom_point()+geom_hline(yintercept=-log10(0.05), col="black",linetype="dotdash")+geom_vline(xintercept=c(-1, 1), col="black",linetype="dotdash")+scale_colour_manual(values = mycolors)+theme_minimal()+theme(axis.text = element_text(size = 12,color = "black"),text =element_text(size = 12,color = "black") )+geom_text_repel() + ggtitle("TcdB1 vs Sham")
ggsave(p, file="./figures/Volcanoplot/TcdB1_vs_Sham_volcano_plot.pdf",height = 5,width = 8)
ggsave(p, file="./figures/Volcanoplot/TcdB1_vs_Sham_volcano_plot.png",height = 5,width = 8)


##### 2 comparison: "TcdB2","Sham"
result2 = results(dds, contrast = c("group","TcdB2","Sham"))
output2= as.data.frame(result2)
# map gene name to gene id
info_tmp=info[rownames(output2),]
output2$gene_name = info_tmp$gene_name
#exporting log2Foldchange table
write.csv(output2, "./log2Foldchange_padj-table/TcdB2-vs-Sham.csv")

# sig gene: abs(log2fc) > 1 && padj<0.05
output2$diffexpressed <- NA
output2$diffexpressed[output2$log2FoldChange > 1 & output2$padj < 0.05] <- "UP"
output2$diffexpressed[output2$log2FoldChange < -1 & output2$padj < 0.05] <- "DOWN"
mycolors <- c("blue", "red", "grey")
names(mycolors) <- c("DOWN", "UP", NA)
output2$delabel <- NA
output2$delabel[!is.na(output2$diffexpressed)] <- output2$gene_name[!is.na(output2$diffexpressed)]
p <- ggplot(data=output2, aes(x=log2FoldChange, y=-log10(padj),color=diffexpressed,label=delabel)) + geom_point()+geom_hline(yintercept=-log10(0.05), col="black",linetype="dotdash")+geom_vline(xintercept=c(-1, 1), col="black",linetype="dotdash")+scale_colour_manual(values = mycolors)+theme_minimal()+theme(axis.text = element_text(size = 12,color = "black"),text =element_text(size = 12,color = "black") )+geom_text_repel() + ggtitle("TcdB2 vs. Sham")
ggsave(p, file="./figures/Volcanoplot/TcdB2_vs_Sham_volcano_plot.pdf",height = 5,width = 8)
ggsave(p, file="./figures/Volcanoplot/TcdB2_vs_Sham_volcano_plot.png",height = 5,width = 8)



##### 3 comparison: "TcdB3","Sham"
result3 = results(dds, contrast = c("group","TcdB3","Sham"))
output3= as.data.frame(result3)
# map gene name to gene id
info_tmp=info[rownames(output3),]
output3$gene_name = info_tmp$gene_name
#exporting log2Foldchange table
write.csv(output3, "./log2Foldchange_padj-table/TcdB3vsSham.csv")

# sig gene: abs(log2fc) > 1 && padj<0.05
output3$diffexpressed <- NA
output3$diffexpressed[output3$log2FoldChange > 1 & output3$padj < 0.05] <- "UP"
output3$diffexpressed[output3$log2FoldChange < -1 & output3$padj < 0.05] <- "DOWN"
mycolors <- c("blue", "red", "grey")
names(mycolors) <- c("DOWN", "UP", NA)
output3$delabel <- NA
output3$delabel[!is.na(output3$diffexpressed)] <- output3$gene_name[!is.na(output3$diffexpressed)]
p <- ggplot(data=output3, aes(x=log2FoldChange, y=-log10(padj),color=diffexpressed,label=delabel)) + geom_point()+geom_hline(yintercept=-log10(0.05), col="black",linetype="dotdash")+geom_vline(xintercept=c(-1, 1), col="black",linetype="dotdash")+scale_colour_manual(values = mycolors)+theme_minimal()+theme(axis.text = element_text(size = 12,color = "black"),text =element_text(size = 12,color = "black") )+geom_text_repel() + ggtitle("TcdB3 vs. Sham")
ggsave(p, file="./figures/Volcanoplot/TcdB3_vs_Sham_volcano_plot.pdf",height = 5,width = 8)
ggsave(p, file="./figures/Volcanoplot/TcdB3_vs_Sham_volcano_plot.png",height = 5,width = 8)




##### 4 comparison: "TcdB1","TcdB3"
result4 = results(dds, contrast = c("group","TcdB1","TcdB3"))
output4= as.data.frame(result4)
# map gene name to gene id
info_tmp=info[rownames(output4),]
output4$gene_name = info_tmp$gene_name
#exporting log2Foldchange table
write.csv(output4, "./log2Foldchange_padj-table/TcdB1vsTcdB3.csv")

# sig gene: abs(log2fc) > 1 && padj<0.05
output4$diffexpressed <- NA
output4$diffexpressed[output4$log2FoldChange > 1 & output4$padj < 0.05] <- "UP"
output4$diffexpressed[output4$log2FoldChange < -1 & output4$padj < 0.05] <- "DOWN"
mycolors <- c("blue", "red", "grey")
names(mycolors) <- c("DOWN", "UP", NA)
output4$delabel <- NA
output4$delabel[!is.na(output4$diffexpressed)] <- output4$gene_name[!is.na(output4$diffexpressed)]
p <- ggplot(data=output4, aes(x=log2FoldChange, y=-log10(padj),color=diffexpressed,label=delabel)) + geom_point()+geom_hline(yintercept=-log10(0.05), col="black",linetype="dotdash")+geom_vline(xintercept=c(-1, 1), col="black",linetype="dotdash")+scale_colour_manual(values = mycolors)+theme_minimal()+theme(axis.text = element_text(size = 12,color = "black"),text =element_text(size = 12,color = "black") )+geom_text_repel() + ggtitle("TcdB1 vs. TcdB3")
ggsave(p, file="./figures/Volcanoplot/TcdB1_vs_TcdB3_volcano_plot.pdf",height = 5,width = 8)
ggsave(p, file="./figures/Volcanoplot/TcdB1_vs_TcdB3_volcano_plot.png",height = 5,width = 8)



##### 5 comparison: "TcdB1","TcdB2"
result4 = results(dds, contrast = c("group","TcdB1","TcdB2"))
output4= as.data.frame(result4)
# map gene name to gene id
info_tmp=info[rownames(output4),]
output4$gene_name = info_tmp$gene_name
#exporting log2Foldchange table
write.csv(output4, "./log2Foldchange_padj-table/TcdB1vsTcdB2.csv")

# sig gene: abs(log2fc) > 1 && padj<0.05
output4$diffexpressed <- NA
output4$diffexpressed[output4$log2FoldChange > 1 & output4$padj < 0.05] <- "UP"
output4$diffexpressed[output4$log2FoldChange < -1 & output4$padj < 0.05] <- "DOWN"
mycolors <- c("blue", "red", "grey")
names(mycolors) <- c("DOWN", "UP", NA)
output4$delabel <- NA
output4$delabel[!is.na(output4$diffexpressed)] <- output4$gene_name[!is.na(output4$diffexpressed)]
p <- ggplot(data=output4, aes(x=log2FoldChange, y=-log10(padj),color=diffexpressed,label=delabel)) + geom_point()+geom_hline(yintercept=-log10(0.05), col="black",linetype="dotdash")+geom_vline(xintercept=c(-1, 1), col="black",linetype="dotdash")+scale_colour_manual(values = mycolors)+theme_minimal()+theme(axis.text = element_text(size = 12,color = "black"),text =element_text(size = 12,color = "black") )+geom_text_repel() + ggtitle("TcdB1 vs. TcdB2")
ggsave(p, file="./figures/Volcanoplot/TcdB1_vs_TcdB2_volcano_plot.pdf",height = 5,width = 8)
ggsave(p, file="./figures/Volcanoplot/TcdB1_vs_TcdB2_volcano_plot.png",height = 5,width = 8)



##### 6 comparison: "TcdB2","TcdB3"
result4 = results(dds, contrast = c("group","TcdB2","TcdB3"))
output4= as.data.frame(result4)
# map gene name to gene id
info_tmp=info[rownames(output4),]
output4$gene_name = info_tmp$gene_name
#exporting log2Foldchange table
write.csv(output4, "./log2Foldchange_padj-table/TcdB2vsTcdB3.csv")

# sig gene: abs(log2fc) > 1 && padj<0.05
output4$diffexpressed <- NA
output4$diffexpressed[output4$log2FoldChange > 1 & output4$padj < 0.05] <- "UP"
output4$diffexpressed[output4$log2FoldChange < -1 & output4$padj < 0.05] <- "DOWN"
mycolors <- c("blue", "red", "grey")
names(mycolors) <- c("DOWN", "UP", NA)
output4$delabel <- NA
output4$delabel[!is.na(output4$diffexpressed)] <- output4$gene_name[!is.na(output4$diffexpressed)]
p <- ggplot(data=output4, aes(x=log2FoldChange, y=-log10(padj),color=diffexpressed,label=delabel)) + geom_point()+geom_hline(yintercept=-log10(0.05), col="black",linetype="dotdash")+geom_vline(xintercept=c(-1, 1), col="black",linetype="dotdash")+scale_colour_manual(values = mycolors)+theme_minimal()+theme(axis.text = element_text(size = 12,color = "black"),text =element_text(size = 12,color = "black") )+geom_text_repel() + ggtitle("TcdB2 vs. TcdB3")
ggsave(p, file="./figures/Volcanoplot/TcdB2_vs_TcdB3_volcano_plot.pdf",height = 5,width = 8)
ggsave(p, file="./figures/Volcanoplot/TcdB2_vs_TcdB3_volcano_plot.png",height = 5,width = 8)



sessionInfo()
```
