require(DESeq2)
require(dplyr)
require(ggplot2)
require(RColorBrewer)
require(pheatmap)
require(viridis)
require(reshape2)

#####
# Analysis of M. tb RNA-seq data from biofilm populations, using DESeq2 and counts from HTSeq.
# Requires: metadata file, HTSeq count files 
#####

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

# subset to only biofilm samples
sampleTable <- sampleTable[sampleTable$Condition == "Biofilm",]

# create sample table for DESeq2 input
DESeq2Table <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = ".",
                                          design = ~ Genotype+Clade)

# remove rows with 0 counts
DESeq2Table <- DESeq2Table[rowSums(counts(DESeq2Table)) > 1,]
DESeq2Table <- DESeq(DESeq2Table)
resultsNames(DESeq2Table)
result <- results(DESeq2Table, alpha = 0.05, lfcThreshold = 0, contrast=c("Genotype","Evolved","Ancestral"))
sum(result$padj <0.05, na.rm = TRUE)
subSet <- subset(result,padj < 0.05) # Supplementary Data 3 - AllPopulations_BFEvA

# BF PCA with arrow
rld <- rlogTransformation(DESeq2Table, blind=FALSE) # log transformation
PCAdata <- DESeq2::plotPCA(rld, intgroup=c("Clade","Genotype"),returnData=T)
write.table(PCAdata,file="DEFiles/BEvA/BF_EvA_PCAdata.tsv",sep="\t")

# PCA data file edited to connect ancestral and evolved samples
PCAdata2 <- read.table("DEFiles/BEvA/BF_EvA_PCAdata_v2.txt",sep="\t",header=T)
rv <- rowVars(assay(rld))
select <- order(rv, decreasing = TRUE)[seq_len(min(500,length(rv)))]
pca <- prcomp(t(assay(rld)[select, ]))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
PC1var <- round(percentVar[1]*100)
PC2var <- round(percentVar[2]*100)

# make plot with arrows (Figure 5A)
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

# separate analysis by lineage 
L4.9.st <- sampleTable[sampleTable$Clade == "L4.9",]
L4.4.st <- sampleTable[sampleTable$Clade == "L4.4",]

# L4.9
L4.9.t <- DESeqDataSetFromHTSeqCount(sampleTable = L4.9.st, directory = ".", design = ~ Genotype)
L4.9.t <- L4.9.t[rowSums(counts(L4.9.t)) > 1,]
L4.9.t <- DESeq(L4.9.t)
L4.9.rld <- rlogTransformation(L4.9.t, blind=FALSE)
DESeq2::plotPCA(L4.9.rld, intgroup=c("Strain","Genotype"))

resultsNames(L4.9.t)
L4.9.r <- results(L4.9.t, alpha = 0.05, lfcThreshold = 0, contrast=c("Genotype","Evolved","Ancestral"))
sum(L4.9.r$padj <0.05, na.rm=TRUE)
L4.9.subset <- subset(L4.9.r, padj<0.05) # Supplementary Data 3 - L4.9_BFEvA

# L4.4
L4.4.t <- DESeqDataSetFromHTSeqCount(sampleTable = L4.4.st, directory = ".", design = ~ Genotype)
L4.4.t <- L4.4.t[rowSums(counts(L4.4.t)) > 1,]
L4.4.t <- DESeq(L4.4.t)
L4.4.rld <- rlogTransformation(L4.4.t, blind=FALSE)
DESeq2::plotPCA(L4.4.rld, intgroup=c("Strain","Genotype"))

resultsNames(L4.4.t)
L4.4.r <- results(L4.4.t, alpha = 0.05, lfcThreshold = 0, contrast=c("Genotype","Evolved","Ancestral"))
sum(L4.4.r$padj <0.05, na.rm=TRUE)
L4.4.subset <- subset(L4.4.r, padj<0.05) # Supplementary Data 3 - L4.4_BFEvA


# volcano plot - L4.9 & L4.4 data plotted together (Figure 5D)
colors <- c(plasma(2,begin=0.4,end = 0.8))

L4.9.voldata <- data.frame(L4.9.r)
L4.4.voldata <- data.frame(L4.4.r)
all.voldata <- rbind(data.frame(padj=L4.9.voldata$padj,lfc=L4.9.voldata$log2FoldChange,gene=rownames(L4.9.voldata),lin="L4.9"),
                     data.frame(padj=L4.4.voldata$padj,lfc=L4.4.voldata$log2FoldChange,gene=rownames(L4.4.voldata),lin="L4.4"))

l4.4.genes <- all.voldata[all.voldata$padj < 0.05 & all.voldata$lin == "L4.4",]$gene
l4.9.genes <- all.voldata[all.voldata$padj < 0.05 & all.voldata$lin == "L4.9",]$gene

all.voldata <- all.voldata %>% mutate(sig=case_when(padj > 0.05 ~ "NS",
                                                    (gene %in% l4.4.genes & gene %in% l4.9.genes) ~ "Both",
                                                    (gene %in% l4.4.genes & !(gene %in% l4.9.genes)) ~ "L4.4",
                                                    (!(gene %in% l4.4.genes) & gene %in% l4.9.genes) ~ "L4.9"))

vol.sub <- all.voldata[all.voldata$sig %in% c("L4.4","L4.9","NS"),]

vol <- ggplot(vol.sub,aes(x=lfc,y=-log10(padj))) + 
  geom_point(aes(color=sig),alpha=0.5,size=2) +
  theme_minimal()+theme(legend.position = "top")+
  scale_color_manual(name="Sub-lineage",labels=c("Both","L4.9","L4.4","NS"),breaks=c("Both","L4.9","L4.4","NS"),values=c("#D76863",colors,"grey"))+
  geom_vline(xintercept=0,linetype="dashed",color="grey27",alpha=0.8) + xlim(-5,5) + ylim(0,8)+
  xlab("log2 Fold Change")

vol

# barplot of up/downregulated genes specific to each lineage (Figure 5C)

shared <- intersect(rownames(L4.9.subset),rownames(L4.4.subset))

c4.4 <- L4.4.subset[rownames(L4.4.subset) %in% shared,]
c4.9 <- L4.9.subset[rownames(L4.9.subset) %in% shared,]
shared.l2fc <- data.frame(gene=rownames(c4.4),l4.4=c4.4$log2FoldChange,l4.9=c4.9$log2FoldChange[match(rownames(c4.4),rownames(c4.9))])
c4.4.u <- L4.4.subset[!(rownames(L4.4.subset) %in% shared),]
c4.9.u <- L4.9.subset[!(rownames(L4.9.subset) %in% shared),]
bar.data <- data.frame(dir=rep(c("up","down"),3), lin=c("L4.4","L4.4","L4.9","L4.9","Shared","Shared"),
                       count=c(sum(c4.4.u$log2FoldChange > 0), sum(c4.4.u$log2FoldChange < 0),
                               sum(c4.9.u$log2FoldChange > 0), sum(c4.9.u$log2FoldChange < 0),
                               sum(shared.l2fc$l4.4 > 0), sum(shared.l2fc$l4.4 < 0)))

count.plot <- ggplot(bar.data,(aes(x=lin,y=count))) + geom_col(aes(fill=lin,alpha=dir),width=0.5)+ theme_bw()+
  scale_fill_manual(name="Sub-lineage",values=c("#FCA636FF","#B12A90FF","#D76863"),labels=c("L4.4","L4.9","Shared"),breaks=c("L4.4","L4.9","Shared"))+
  scale_alpha_manual(name="DEGs",values=c(0.4,1),labels=c("downregulated","upregulated"),breaks=c("down","up"))+
  xlab(NULL)+ylab("Count")+
  guides(fill=guide_legend(order=1))+
  theme(axis.text=element_text(size=12),axis.title.y = element_text(size=12))
count.plot

# separate by individual strains (Supplementary Data 3 - BFEvA for each strain)
strains <- c("31","49","55","72","345","540")

for (strain in strains){
  st <- sampleTable[sampleTable$Strain == strain,]
  t <- DESeqDataSetFromHTSeqCount(sampleTable = st, directory = ".", design = ~ Genotype)
  t <- t[rowSums(counts(t)) > 1,]
  t <- DESeq(t)
  r <- results(t, alpha = 0.05, lfcThreshold = 0, contrast=c("Genotype","Evolved","Ancestral"))
  sub <- subset(r, padj < 0.05)
  write.csv(sub,paste("../2023.02.28_newFigs/DEFiles/",paste(strain,"BFEvA",sep="_"),".csv",sep=""))
}

# l2fc values by feature type (Figure 8B)
ncrna <- read.delim("Metadata/all_ncRNA.txt",header=F)$V1
sorfs <- read.delim("Metadata/all_sORFs.txt",header=F)$V1
result <- data.frame(result)
result$gene <- rownames(result)

result <- result %>% mutate(type=case_when(gene %in% ncrna ~ "ncRNA",
                                           gene %in% sorfs ~ "sORF",
                                           TRUE ~ "ORF"))
my_colors <- magma(n=4, end = 0.8, begin=0.25)
feature_colors <- c(my_colors[c(2)],"forestgreen",my_colors[c(4)])

plot <- ggplot(result,aes(x=type,y=log2FoldChange)) + geom_hline(yintercept=0,linetype="dashed")+
  geom_boxplot(aes(fill=type),alpha=0.8)+
  theme_minimal() + xlab(NULL)+ ylab("Log2 Fold Change")+
  scale_fill_manual(name=NULL,values=feature_colors,breaks=c("ncRNA","ORF","sORF"))+
  theme(axis.text=element_text(size=12,face="bold"),
        axis.text.x = element_blank(),
        axis.title = element_text(size=12),
        strip.text = element_text(size=12,face="bold"))
plot

stat <- compare_means(log2FoldChange ~ type, data = result,p.adjust.method = "BH")
stat <- stat %>% mutate(y.position=c(5,NA,6))
box.stats <- plot + stat_pvalue_manual(stat, label="p.signif",label.size=4)
box.stats

ExportPlot(box.stats,"../2023.02.28_newFigs/IndvFigures/BEvA_featuretype",width=6,height=6)
