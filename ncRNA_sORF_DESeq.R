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

setwd("~/Desktop/2023.04.26_RNAseq/mtb_ExpEvo_RNA/")

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
E.table <- sampleTable[sampleTable$Genotype == "Evolved",]

# DE calculations by condition (biofilm vs plantkonic)
E.DESeq <- DESeqDataSetFromHTSeqCount(sampleTable = E.table, directory = ".", design = ~ Condition)

# filter out genes with 0 counts in all samples
E.DESeq <- E.DESeq[rowSums(counts(E.DESeq)) > 1,]
nrow(E.DESeq)

# only include ncRNA and sORFs
ncRNA.sORF.DESeq <- E.DESeq[rownames(E.DESeq) %in% ncrna_sorfs$feature,]
nrow(ncRNA.sORF.DESeq)

# run DESeq analysis
ncRNA.sORF.DESeq <- DESeq(ncRNA.sORF.DESeq)
resultsNames(ncRNA.sORF.DESeq)

# r-log transformation
ncRNA.sORF.rld <- rlogTransformation(ncRNA.sORF.DESeq, blind=FALSE)

# results of biofilm vs planktonic comparison
ncRNA.sORF.result <- data.frame(results(ncRNA.sORF.DESeq, alpha = 0.05,lfcThreshold = 0, 
                                  contrast=c("Condition","Biofilm","Planktonic")))
ncRNA.sORF.result$gene <- rownames(ncRNA.sORF.result)

ncRNA.sORF.result <- ncRNA.sORF.result %>% mutate(type=case_when(gene %in% ncrna ~ "ncRNA",
                                           gene %in% sorfs ~ "sORF"))

# volcano plot - ncRNA & sORF (Figure 7A)
voldata <- ncRNA.sORF.result
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
        axis.title = element_text(size=12), strip.text = element_text(size=12, face="bold"),
        legend.text = element_text(size=12)) + 
  facet_wrap(~type,ncol=1) + xlab("log2 fold change")

vol

ExportPlot(vol,"../NewFigures/Figure7A",width=6,height=8)


# heatmap of sDE ncRNA and sORFs (Figure 7B)
ncRNA.sORF.result.sig <- ncRNA.sORF.result[ncRNA.sORF.result$padj < 0.05,]
BvPCounts <- assay(ncRNA.sORF.rld)[ncRNA.sORF.result.sig$gene,]
BvPCounts <- BvPCounts - rowMeans(BvPCounts)

my_colors <- magma(n=4, end = 0.8, begin=0.1)
condition <- my_colors[c(1,3)]
type <- my_colors[c(2,4)]
custom.palette <- c("red","white","blue")
colors <- colorRampPalette(rev(custom.palette))(20)
names(condition) <- c("Biofilm","Planktonic")
names(type) <- c("ncRNA","sORF")
anno_colors <- list(Condition=condition,Type=type)
anno.c <- as.data.frame(colData(ncRNA.sORF.rld)[c("Condition")])
anno.r <- data.frame("Type"= ncRNA.sORF.result.sig$type)
rownames(anno.r) <- ncRNA.sORF.result.sig$gene

de <- pheatmap.type(BvPCounts, annRow = anno.r,color = colors, show_rownames = F, annotation_col = anno.c,
              show_colnames = F, cutree_cols = 2, cutree_rows=2, annotation_colors = anno_colors,
              kmeans_k = NA, border_color = NA, scale= "column", fontsize = 12)


# PCA of evolved and ancestral biofilm samples - ncRNA & sORF only (Figure 7C)

# subset to only biofilm samples
B.table <- sampleTable[sampleTable$Condition == "Biofilm",]

# create sample table for DESeq2 input
B.DESeq <- DESeqDataSetFromHTSeqCount(sampleTable = B.table, directory = ".",
                                          design = ~ Genotype+Clade)

# remove rows with 0 counts
B.DESeq <- B.DESeq[rowSums(counts(B.DESeq)) > 1,]

# only include ncRNA and sORFs
B.DESeq <- B.DESeq[rownames(B.DESeq) %in% ncrna_sorfs$feature,]

# log transformation
rld <- rlogTransformation(B.DESeq, blind=FALSE)

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
  theme(legend.position = "top",axis.text=element_text(size=12),axis.title=element_text(size=12),
        legend.text=element_text(size=12))

PCAplot

ExportPlot(PCAplot,"../NewFigures/Figure7C",width=12,height=6)

# PCA of evolved and ancestral planktonic samples - ncRNA & sORF only (Figure S6)


# subset to only biofilm samples
P.table <- sampleTable[sampleTable$Condition == "Planktonic",]

# create sample table for DESeq2 input
DESeq2Table <- DESeqDataSetFromHTSeqCount(sampleTable = P.table, directory = ".",
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
  theme(legend.position = "top",axis.text=element_text(size=12),axis.title=element_text(size=12),
        legend.text=element_text(size=12))

PCAplot

ExportPlot(PCAplot,"../NewFigures/Supplement/FigureS6",width=12,height=6)

