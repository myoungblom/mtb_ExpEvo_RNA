require(DESeq2)
require(dplyr)
require(ggplot2)
require(RColorBrewer)
require(pheatmap)
require(viridis)
require(reshape2)
require(ggpubr)
require(Pigengene)

#####
# Analysis of ncRNA and sORF expression from evolved populations, using DESeq2 and counts from HTSeq.
# Requires: metadata file, HTSeq count files 
#####

# lists of ncRNA and sORFs
ncrna <- read.delim("Metadata/all_ncRNA.txt",header=F)$V1
sorfs <- read.delim("Metadata/all_sORFs.txt",header=F)$V1
ncrna_sorfs <- rbind(data.frame(feature=ncrna,type="ncRNA"),
                     data.frame(feature=sorfs,type="sORF"))

# read metadata file
sampleData <- read.delim("Metadata/mtb_ExpEvo_RNA_metadata.txt", sep="", header = TRUE)

# reformat metadata
sampleTable <- data.frame(SampleName = sampleData$LibraryName,
                          CountsFile = sampleData$CountsFile,
                          Strain = factor(sampleData$Strain),
                          Genotype = factor(sampleData$Genotype),
                          Condition = factor(sampleData$Condition,levels = c("Planktonic","Biofilm")),
                          WetWeight = as.numeric(sampleData$WetWeight),
                          Clade = factor(sampleData$Clade),
                          SampleID = sampleData$SampleID)

# subset to only evolved samples
sampleTable <- sampleTable[sampleTable$Genotype == "Evolved",]

# DE calculations by condition (biofilm vs plantkonic)
DESeq2Table <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = ".", design = ~ Condition)

# filter out genes with 0 counts in all samples
DESeq2Table <- DESeq2Table[rowSums(counts(DESeq2Table)) > 1,]
nrow(DESeq2Table)

# only include ncRNA and sORFs
DESeq2Table <- DESeq2Table[rownames(DESeq2Table) %in% ncrna_sorfs$feature,]
nrow(DESeq2Table)

# run DESeq analysis
DESeq2Table <- DESeq(DESeq2Table)
resultsNames(DESeq2Table)

# r-log transformation
rld <- rlogTransformation(DESeq2Table, blind=FALSE)

# results of biofilm vs planktonic comparison
result <- data.frame(results(DESeq2Table, alpha = 0.05,lfcThreshold = 0, contrast=c("Condition","Biofilm","Planktonic")))
result$gene <- rownames(result)

result <- result %>% mutate(type=case_when(gene %in% ncrna ~ "ncRNA",
                                           gene %in% sorfs ~ "sORF"))

# volcano plot - ncRNA & sORF (Figure 7A)
voldata <- result
voldata <- voldata %>% mutate(sig=case_when(padj < 0.05 & log2FoldChange > 0 ~ "Significant-Up",
                                            padj < 0.05 & log2FoldChange < 0 ~ "Significant-Down",
                                            padj > 0.05 ~ "Not significant"))


custom.palette <- colorRampPalette(rev(c("red","white","blue")))(20)
my_colors <- c("grey",custom.palette[15],custom.palette[5])

vol <- ggplot(voldata,aes(x=log2FoldChange,y=-log10(padj))) + geom_point(aes(color=sig),alpha=0.8) + 
  theme_minimal()+xlim(-11,11)+scale_color_manual(name=NULL,values=my_colors,
                                                  breaks=c("Not significant","Significant-Up","Significant-Down"),
                                                  labels=c("Not significant","Upregulated","Downregulated"))+
  theme(legend.position = "top",axis.text = element_text(size=10),
        axis.title = element_text(size=10), strip.text = element_text(size=12, face="bold")) + 
  facet_wrap(~type,ncol=1) + xlab("log2 fold change")

vol

ExportPlot(vol,"../2023.02.28_newFigs/IndvFigures/ncRNA_sORFs_volcano",width=6,height=8)


# heatmap of sDE ncRNA and sORFs (Figure 7B)
result.sig <- result[result$padj < 0.05,]
BvPCounts <- assay(rld)[result.sig$gene,]
BvPCounts <- BvPCounts - rowMeans(BvPCounts)

my_colors <- magma(n=4, end = 0.8, begin=0.1)
condition <- my_colors[c(1,3)]
type <- my_colors[c(2,4)]
custom.palette <- c("red","white","blue")
colors <- colorRampPalette(rev(custom.palette))(20)
names(condition) <- c("Biofilm","Planktonic")
names(type) <- c("ncRNA","sORF")
anno_colors <- list(Condition=condition,Type=type)
anno.c <- as.data.frame(colData(rld)[c("Condition")])
anno.r <- data.frame("Type"= result.sig$type)
rownames(anno.r) <- result.sig$gene

de <- pheatmap.type(BvPCounts, annRow = anno.r,color = colors, show_rownames = F, annotation_col = anno.c,
              show_colnames = F, cutree_cols = 2, cutree_rows=2, annotation_colors = anno_colors,
              kmeans_k = NA, border_color = NA, scale= "column", fontsize = 12)


ExportPlot(de.heatmap,"../2023.02.28_newFigs/IndvFigures/ncRNA_sORFs_heatmap",width=8,height=10)

# PCA of evolved and ancestral biofilm samples - ncRNA & sORF only (Figure 7C)

# reformat metadata
sampleTable <- data.frame(SampleName = sampleData$LibraryName,
                          CountsFile = sampleData$CountsFile,
                          Strain = factor(sampleData$Strain),
                          Genotype = factor(sampleData$Genotype),
                          Condition = factor(sampleData$Condition,levels = c("Planktonic","Biofilm")),
                          WetWeight = as.numeric(sampleData$WetWeight),
                          Clade = factor(sampleData$Clade),
                          SampleID = sampleData$SampleID)

# subset to only biofilm samples
sampleTable <- sampleTable[sampleTable$Condition == "Biofilm",]

# create sample table for DESeq2 input
DESeq2Table <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = ".",
                                          design = ~ Genotype+Clade)

# remove rows with 0 counts
DESeq2Table <- DESeq2Table[rowSums(counts(DESeq2Table)) > 1,]

# only include ncRNA and sORFs
DESeq2Table <- DESeq2Table[rownames(DESeq2Table) %in% ncrna_sorfs$feature,]

# log transformation
rld <- rlogTransformation(DESeq2Table, blind=FALSE)

PCAdata <- DESeq2::plotPCA(rld, intgroup=c("Clade","Genotype"),returnData=T)
write.table(PCAdata,file="DEFiles/BEvA/BF_EvA_PCAdata_ncRNAsORFs.tsv",sep="\t")

# PCA data file edited to connect ancestral and evolved samples
PCAdata2 <- read.table("DEFiles/BEvA/BF_EvA_PCAdata_ncRNAsORFs_v2.txt",sep="\t",header=T)
rv <- rowVars(assay(rld))
select <- order(rv, decreasing = TRUE)[seq_len(min(500,length(rv)))]
pca <- prcomp(t(assay(rld)[select, ]))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
PC1var <- round(percentVar[1]*100)
PC2var <- round(percentVar[2]*100)

colors <- c(plasma(2,begin=0.4,end = 0.8))
arrows <- arrow(ends="last", length=unit(0.1,"inches"),angle=20,type="closed")

PCAplot <- ggplot()+geom_point(data=PCAdata2,aes(x=PC1,y=PC2,color=Clade,shape=Genotype),size=4,alpha=0.9)+
  geom_curve(aes(x=PC1,y=PC2,xend=PC12,yend=PC22,color=Clade),alpha=0.9,data=PCAdata2,
             arrow=arrows,curvature = 0.3,na.rm = T,show.legend = F)+
  ylim(-25,25)+theme_minimal()+xlab(paste0("PC1: ",PC1var,"%"))+ylab(paste0("PC2: ",PC2var,"%"))+
  scale_color_manual(name="Sub-lineage",labels=c("L4.9","L4.4"),breaks=c("L4.9","L4.4"),values=colors)+
  scale_shape_manual(name="Genotype",breaks=c("Ancestral","Evolved"), values=c(10,16))+
  theme(legend.position = "top",axis.text=element_text(size=8),axis.title=element_text(size=10),
        legend.text=element_text(size=12))

PCAplot

ExportPlot(PCAplot,"../2023.02.28_newFigs/IndvFigures/ncRNA_sORFs_BEvA_PCA",width=14,height=8)

# PCA of evolved and ancestral planktonic samples - ncRNA & sORF only (Figure S6)

# reformat metadata
sampleTable <- data.frame(SampleName = sampleData$LibraryName,
                          CountsFile = sampleData$CountsFile,
                          Strain = factor(sampleData$Strain),
                          Genotype = factor(sampleData$Genotype),
                          Condition = factor(sampleData$Condition,levels = c("Planktonic","Biofilm")),
                          WetWeight = as.numeric(sampleData$WetWeight),
                          Clade = factor(sampleData$Clade),
                          SampleID = sampleData$SampleID)

# subset to only biofilm samples
sampleTable <- sampleTable[sampleTable$Condition == "Planktonic",]

# create sample table for DESeq2 input
DESeq2Table <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = ".",
                                          design = ~ Genotype+Clade)

# remove rows with 0 counts
DESeq2Table <- DESeq2Table[rowSums(counts(DESeq2Table)) > 1,]

# only include ncRNA and sORFs
DESeq2Table <- DESeq2Table[rownames(DESeq2Table) %in% ncrna_sorfs$feature,]

# log transformation
rld <- rlogTransformation(DESeq2Table, blind=FALSE)

PCAdata <- DESeq2::plotPCA(rld, intgroup=c("Clade","Genotype"),returnData=T)
write.table(PCAdata,file="DEFiles/PEvA/P_EvA_PCAdata_ncRNAsORFs.tsv",sep="\t")

# PCA data file edited to connect ancestral and evolved samples
PCAdata2 <- read.table("DEFiles/PEvA/P_EvA_PCAdata_ncRNAsORFs_v2.txt",sep="\t",header=T)
rv <- rowVars(assay(rld))
select <- order(rv, decreasing = TRUE)[seq_len(min(500,length(rv)))]
pca <- prcomp(t(assay(rld)[select, ]))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
PC1var <- round(percentVar[1]*100)
PC2var <- round(percentVar[2]*100)

colors <- c(plasma(2,begin=0.4,end = 0.8))
arrows <- arrow(ends="last", length=unit(0.1,"inches"),angle=20,type="closed")

PCAplot <- ggplot()+geom_point(data=PCAdata2,aes(x=PC1,y=PC2,color=Clade,shape=Genotype),size=4,alpha=0.9)+
  geom_curve(aes(x=PC1,y=PC2,xend=PC12,yend=PC22,color=Clade),alpha=0.9,data=PCAdata2,
             arrow=arrows,curvature = 0.3,na.rm = T,show.legend = F)+
  ylim(-25,25)+theme_minimal()+xlab(paste0("PC1: ",PC1var,"%"))+ylab(paste0("PC2: ",PC2var,"%"))+
  scale_color_manual(name="Sub-lineage",labels=c("L4.9","L4.4"),breaks=c("L4.9","L4.4"),values=colors)+
  scale_shape_manual(name="Genotype",breaks=c("Ancestral","Evolved"), values=c(10,16))+
  theme(legend.position = "top",axis.text=element_text(size=8),axis.title=element_text(size=10),
        legend.text=element_text(size=12))

PCAplot

ExportPlot(PCAplot,"../2023.02.28_newFigs/IndvFigures/ncRNA_sORFs_PEvA_PCA",width=14,height=8)

# boxplot
result.sig <- result[result$padj < 0.05,]
result.sig$type <- factor(result.sig$type,levels=c("ncRNA","sORF","ORF"))
box <- ggplot(result.sig,aes(x=type,y=log2FoldChange)) + geom_hline(yintercept=0,linetype="dashed")+
  geom_boxplot(aes(fill=type)) + theme_minimal()+
  scale_fill_manual(name=NULL,values=c(my_colors[c(2,4)],"#339966"),breaks=c("ncRNA","sORF","ORF"))+
  xlab(NULL) + theme(legend.position = "none",axis.text.x=element_text(size=12,face="bold"))

stat <- compare_means(log2FoldChange ~ type, data = result.sig, p.adjust.method = "BH")
stat <- stat %>% mutate(y.position=c(4,5,6))
box.stats <- box + stat_pvalue_manual(stat, label="p.signif",label.size=4)
box.stats

ExportPlot(box.stats,"Figures/feature_l2fc_boxplots",width=4,height=4)

e.result <- result
e.result <- e.result %>% mutate(type=case_when(X %in% ncrna ~ "ncRNA", X %in% sorfs ~ "sORF",TRUE ~ "ORF"))

a.result <- as.data.frame(read.csv("../DE_AncBvP/DElists/Anc_BvP_allGenes.csv"))
a.result <- a.result %>% mutate(type=case_when(X %in% ncrna ~ "ncRNA", X %in% sorfs ~ "sORF",TRUE ~ "ORF"))

all <- rbind(data.frame(genotype="Ancestral",l2fc=a.result$log2FoldChange,p=a.result$padj,type=a.result$type),
             data.frame(genotype="Evolved",l2fc=e.result$log2FoldChange,p=e.result$padj,type=e.result$type))

all.sig <- subset(all,p<0.05)
box <- ggplot(all,aes(x=type,y=l2fc)) + geom_hline(yintercept=0,linetype="dashed")+
  geom_boxplot(aes(fill=type)) + theme_minimal()+
  scale_fill_manual(name=NULL,values=c(my_colors[c(2,4)],"#339966"),breaks=c("ncRNA","sORF","ORF"))+
  xlab(NULL) + theme(legend.position = "none",axis.text.x=element_text(size=12,face="bold"))+
  facet_wrap(~genotype)

stat <- compare_means(l2fc ~ type, data = all, group.by = "genotype",p.adjust.method = "BH")
stat <- stat %>% mutate(y.position=c(4.5,5.5,6.5,4.5,5.5,6.5))
box.stats <- box + stat_pvalue_manual(stat, label="p.signif",label.size=4)
box.stats


#  Evolved vs Ancestral
setwd("~/Desktop/2022.12.15_sORF_RNAseq/")
sampleData <- read.delim("./metadata/2022_allRNA_metadata.txt", sep="", header = TRUE)

sampleTable <- data.frame(SampleName = sampleData$LibraryName,
                          CountsFile = sampleData$CountsFile,
                          Type = factor(sampleData$Type),
                          Strain = factor(sampleData$Strain),
                          Genotype = factor(sampleData$Genotype),
                          Condition = factor(sampleData$Condition,levels = c("Planktonic","Biofilm")),
                          Clade = factor(sampleData$Clade),
                          SampleID = sampleData$SampleID)

sampleTable <- sampleTable[sampleTable$Type == "G" & sampleTable$Condition == "Planktonic",]

DESeq2Table <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = ".", design = ~ Genotype)
DESeq2Table <- DESeq2Table[rowSums(counts(DESeq2Table)) > 1,]
ncRNATable <- DESeq2Table[rownames(DESeq2Table) %in% ncrna | rownames(DESeq2Table) %in% sorfs]

ncRNATable <- DESeq(ncRNATable)

b.rld <- rlogTransformation(ncRNATable, blind=FALSE)
p.rld <- rlogTransformation(ncRNATable, blind=FALSE)

b.results <- results(ncRNATable, alpha = 0.05, lfcThreshold = 1, contrast=c("Genotype","Evolved","Ancestral"))
b.results$gene <- rownames(b.results)
p.results <- results(ncRNATable, alpha = 0.05, lfcThreshold = 1, contrast=c("Genotype","Evolved","Ancestral"))
p.results$gene <- rownames(p.results)

ncrna <- read.delim("DE_ncRNA_sORF/all_ncRNA.txt",header=F)$V1
sorfs <- read.delim("DE_ncRNA_sORF/all_sORFs.txt",header=F)$V1

bp.results <- rbind(data.frame(l2fc=b.results$log2FoldChange,condition="Biofilm",gene=b.results$gene,p=b.results$padj),
                    data.frame(l2fc=p.results$log2FoldChange,condition="Planktonic",gene=p.results$gene,p=p.results$padj))

bp.results <- bp.results %>% mutate(type=case_when(gene %in% ncrna ~ "ncRNA",
                                           gene %in% sorfs ~ "sORF",
                                           TRUE ~ "ORF"))
bp.results.sig <- subset(bp.results,p<0.05)
bp.results.sig$type <- factor(bp.results.sig$type,levels=c("ncRNA","sORF","ORF"))

my_colors <- magma(n=4, end = 0.8, begin=0.1)

box <- ggplot(bp.results,aes(x=type,y=l2fc)) + geom_hline(yintercept=0,linetype="dashed")+
  geom_boxplot(aes(fill=type)) + theme_minimal()+
  scale_fill_manual(name=NULL,values=c(my_colors[c(2,4)],"#339966"),breaks=c("ncRNA","sORF","ORF"))+
  xlab(NULL) + theme(legend.position = "none",axis.text.x=element_text(size=12,face="bold"))+
  facet_wrap(~condition)+ylab("log2 fold change (evolved vs ancestral)")
box

stat <- compare_means(l2fc ~ type, data = bp.results.sig, group.by = "condition",p.adjust.method = "BH")
stat <- stat %>% mutate(y.position=c(NA,3,NA))
box.stats <- box + stat_pvalue_manual(stat, label="p.adj",label.size=4)
box.stats

ExportPlot(box.stats,"DE_ncRNA_sORF/Figures/feature_l2fc_boxplots_EvA",width=6,height=4)

# PCA with arrows
#b.rld.sub <- b.rld[rownames(b.rld) %in% ncrna | rownames(b.rld) %in% sorfs,]
p.rld.sub <- p.rld[rownames(p.rld) %in% ncrna | rownames(p.rld) %in% sorfs,]

DESeq2::plotPCA(b.rld, intgroup=c("Clade","Genotype","Condition"))

PCAdata <- DESeq2::plotPCA(b.rld.sub, intgroup=c("Clade","Genotype"),returnData=T)
#write.table(PCAdata,file="DE_ncRNA_sORF/BFEvA_ncRNA_PCAdata_v3.tsv",sep="\t")

PCAdata2 <- read.table("DE_ncRNA_sORF/BFEvA_ncRNA_PCAdata_v4.txt",sep="\t",header=T)
rv <- rowVars(assay(b.rld.sub))
select <- order(rv, decreasing = TRUE)[seq_len(min(500,length(rv)))]
pca <- prcomp(t(assay(b.rld.sub)[select, ]))
percentVar <- pca$sdev^2/sum(pca$sdev^2);percentVar
PC1var <- round(percentVar[1]*100)
PC2var <- round(percentVar[2]*100)

colors <- c(plasma(2,begin=0.4,end = 0.8))
arrows <- arrow(ends="last", length=unit(0.1,"inches"),angle=20,type="closed")

PCAplot <- ggplot()+geom_point(data=PCAdata2,aes(x=PC1,y=PC2,color=Clade,shape=Genotype),size=4,alpha=0.9)+
  geom_curve(aes(x=PC1,y=PC2,xend=PC12,yend=PC22,color=Clade),alpha=0.9,data=PCAdata2,
             arrow=arrows,curvature = 0.3,na.rm = T,show.legend = F)+
  ylim(-25,25)+theme_minimal()+xlab(paste0("PC1: ",PC1var,"%"))+ylab(paste0("PC2: ",PC2var,"%"))+
  scale_color_manual(name="Sub-Lineage",labels=c("L4.9","L4.4"),breaks=c("C","JM_DS6"),values=colors)+
  scale_shape_manual(name="Genotype",breaks=c("Ancestral","Evolved"), values=c(10,16))+
  theme(legend.position = "top",axis.text=element_text(size=8),axis.title=element_text(size=10),
        legend.text=element_text(size=12))

PCAplot

ExportPlot(PCAplot,"DE_ncRNA_sORF/Figures/ncRNA_BFEvA_PCAwArrows_v3",width=14,height=8)

ExportPlot <- function(gplot, filename, width=2, height=1.5) {
  # Export plot in PDF and EPS.
  # Notice that A4: width=11.69, height=8.27
  ggsave(paste(filename, '.pdf', sep=""), gplot, width = width, height = height)
  postscript(file = paste(filename, '.eps', sep=""), width = width, height = height, family = "sans")
  print(gplot)
  dev.off()
  png(file = paste(filename, '_.png', sep=""), width = width * 100, height = height * 100)
  print(gplot)
  dev.off()
}
