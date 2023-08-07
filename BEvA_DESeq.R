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

#function to export plot
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

setwd("~/Desktop/2023.05.02_RNAseq/mtb_ExpEvo_RNA/")

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
B.table <- sampleTable[sampleTable$Condition == "Biofilm",]

# create sample table for DESeq2 input
B.DESeq <- DESeqDataSetFromHTSeqCount(sampleTable = B.table, directory = ".",design = ~ Genotype)

# remove rows with 0 counts
B.DESeq <- B.DESeq[rowSums(counts(B.DESeq)) > 1,]
B.DESeq <- DESeq(B.DESeq)
resultsNames(B.DESeq)
BEvA.result <- results(B.DESeq, alpha = 0.05, lfcThreshold = 0, contrast=c("Genotype","Evolved","Ancestral"))
#write.csv(x = BEvA.result,"DEFiles/BEvA/all_BEvA.csv",quote=F)

sum(BEvA.result$padj <0.05, na.rm = TRUE)
subSet <- subset(BEvA.result,padj < 0.05) # Supplementary Data 3 - AllPopulations_BFEvA
#write.csv(x = subSet,"DEFiles/BEvA/all_BEvA_DE.csv",quote=F)

# BF PCA with arrow
B.rld <- rlogTransformation(B.DESeq, blind=FALSE) # log transformation
PCAdata <- DESeq2::plotPCA(B.rld, intgroup=c("Clade","Genotype"),returnData=T)
#write.table(PCAdata,file="DEFiles/BEvA/BF_EvA_PCAdata.tsv",sep="\t")

# PCA data file edited to connect ancestral and evolved samples
PCAdata2 <- read.table("DEFiles/BEvA/BF_EvA_PCAdata_v2.txt",sep="\t",header=T)
rv <- rowVars(assay(B.rld))
select <- order(rv, decreasing = TRUE)[seq_len(min(500,length(rv)))]
pca <- prcomp(t(assay(B.rld)[select, ]))
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
  theme(legend.position = "top",axis.text=element_text(size=10),axis.title=element_text(size=12),
        legend.text=element_text(size=12))

PCAplot

ExportPlot(PCAplot,"../NewFigures/Figure5A",width=6,height=5)

# separate analysis by lineage 
L4.9.st <- B.table[B.table$Clade == "L4.9",]
L4.4.st <- B.table[B.table$Clade == "L4.4",]

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
write.csv(L4.9.subset,"DEFiles/BEvA/L4.9_BEvA_DE.csv",quote=F)

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
write.csv(L4.4.subset,"DEFiles/BEvA/L4.4_BEvA_DE.csv",quote=F)

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
  geom_point(aes(color=sig),alpha=0.75,size=2) +
  theme_minimal()+theme(legend.position = "top")+
  scale_color_manual(name="Sub-lineage",labels=c("Both","L4.9","L4.4","NS"),breaks=c("Both","L4.9","L4.4","NS"),values=c("#D76863",colors,"grey"))+
  geom_vline(xintercept=0,linetype="dashed",color="grey27",alpha=0.8) + xlim(-5,5) + ylim(0,8)+
  xlab("log2 Fold Change") +
  theme(legend.position = "top",axis.text=element_text(size=10),axis.title=element_text(size=12),
        legend.text=element_text(size=12))

vol

ExportPlot(vol,"../NewFigures/Figure5D",width=6,height=6)

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
  scale_alpha_manual(name=NULL,values=c(0.4,1),labels=c("downregulated","upregulated"),breaks=c("down","up"))+
  xlab(NULL)+ylab("Count")+
  guides(fill=guide_legend(order=1))+
  theme(axis.text=element_text(size=12),axis.title.y = element_text(size=12))
count.plot

ExportPlot(count.plot,"../NewFigures/Figure5C",width=6,height=6)

l4.9.up <- (nrow(subset(c4.9.u,log2FoldChange > 0))/nrow(c4.9.u))*100
l4.4.up <- (nrow(subset(c4.4.u,log2FoldChange > 0))/nrow(c4.4.u))*100

print(paste("% of unique L4.9 genes that are upregulated: ",l4.9.up,"%"))
print(paste("% of unique L4.4 genes that are upregulated: ",l4.4.up,"%"))

# separate by individual strains (Supplementary Data 3 - BFEvA for each strain)
strains <- c("31","49","55","72","345","540")

for (strain in strains){
  st <- B.table[B.table$Strain == strain,]
  t <- DESeqDataSetFromHTSeqCount(sampleTable = st, directory = ".", design = ~ Genotype)
  t <- t[rowSums(counts(t)) > 1,]
  t <- DESeq(t)
  r <- results(t, alpha = 0.05, lfcThreshold = 0, contrast=c("Genotype","Evolved","Ancestral"))
  write.csv(r,paste("DEFiles/BEvA/",paste(strain,"BFEvA","allGenes",sep="_"),".csv",sep=""),quote=F)
  sub <- subset(r, padj < 0.05)
  write.csv(sub,paste("DEFiles/BEvA/",paste(strain,"BFEvA",sep="_"),".csv",sep=""),quote=F)
}

# l2fc values by feature type (Figure 8B)
ncrna <- read.delim("Metadata/all_ncRNA.txt",header=F)$V1
sorfs <- read.delim("Metadata/all_sORFs.txt",header=F)$V1
BEvA.result <- data.frame(BEvA.result)
BEvA.result$gene <- rownames(BEvA.result)

BEvA.result <- BEvA.result %>% mutate(type=case_when(gene %in% ncrna ~ "ncRNA",
                                           gene %in% sorfs ~ "sORF",
                                           TRUE ~ "ORF"))
my_colors <- magma(n=4, end = 0.8, begin=0.25)
feature_colors <- c(my_colors[c(2)],"forestgreen",my_colors[c(4)])

plot <- ggplot(BEvA.result,aes(x=type,y=log2FoldChange)) + geom_hline(yintercept=0,linetype="dashed")+
  geom_boxplot(aes(fill=type),alpha=0.8)+
  theme_minimal() + xlab(NULL)+ ylab("Log2 Fold Change")+
  scale_fill_manual(name=NULL,values=feature_colors,breaks=c("ncRNA","ORF","sORF"))+
  theme(axis.text=element_text(size=12,face="bold"),
        axis.text.x = element_blank(),
        axis.title = element_text(size=12),
        strip.text = element_text(size=12,face="bold"))+
  scale_y_continuous(limits=c(-13,6),breaks=c(-10,-5,0,5))
plot

stat <- compare_means(log2FoldChange ~ type, data = BEvA.result,p.adjust.method = "BH")
stat <- stat %>% mutate(y.position=c(5,NA,6))
box.stats <- plot + stat_pvalue_manual(stat, label="p.signif",label.size=4)
box.stats

ExportPlot(box.stats,"../NewFigures/Figure8B",width=6,height=6)

# lpdA operon expression (Figure S6)
operon.genes <- c("lpdA","glpD2","phoY1","Rv3300c","atsB","lpqC")
b.540 <- read.csv("DEFiles/BEvA/540_BFEvA_allGenes.csv")
p.540 <- read.csv("DEFiles/PEvA/540_PEvA_allGenes.csv")
b.49 <- read.csv("DEFiles/BEvA/49_BFEvA_allGenes.csv")
p.49 <- read.csv("DEFiles/PEvA/49_PEvA_allGenes.csv")

b.345 <- read.csv("DEFiles/BEvA/345_BFEvA_allGenes.csv")
p.345 <- read.csv("DEFiles/PEvA/345_PEvA_allGenes.csv")
b.72 <- read.csv("DEFiles/BEvA/72_BFEvA_allGenes.csv")
p.72 <- read.csv("DEFiles/PEvA/72_PEvA_allGenes.csv")

b.31 <- read.csv("DEFiles/BEvA/31_BFEvA_allGenes.csv")
p.31 <- read.csv("DEFiles/PEvA/31_PEvA_allGenes.csv")
b.55 <- read.csv("DEFiles/BEvA/55_BFEvA_allGenes.csv")
p.55 <- read.csv("DEFiles/PEvA/55_PEvA_allGenes.csv")

all.lpda <- rbind(data.frame(gene=b.540$X,l2fc=b.540$log2FoldChange,pvalue=b.540$padj,condition="Biofilm",strain="MT540"),
                  data.frame(gene=p.540$X,l2fc=p.540$log2FoldChange,pvalue=p.540$padj,condition="Planktonic",strain="MT540"),
                  data.frame(gene=b.49$X,l2fc=b.49$log2FoldChange,pvalue=b.49$padj,condition="Biofilm",strain="MT49"),
                  data.frame(gene=p.49$X,l2fc=p.49$log2FoldChange,pvalue=p.49$padj,condition="Planktonic",strain="MT49"),
                  data.frame(gene=b.345$X,l2fc=b.345$log2FoldChange,pvalue=b.345$padj,condition="Biofilm",strain="MT345"),
                  data.frame(gene=p.345$X,l2fc=p.345$log2FoldChange,pvalue=p.345$padj,condition="Planktonic",strain="MT345"),
                  data.frame(gene=b.72$X,l2fc=b.72$log2FoldChange,pvalue=b.72$padj,condition="Biofilm",strain="MT72"),
                  data.frame(gene=p.72$X,l2fc=p.72$log2FoldChange,pvalue=p.72$padj,condition="Planktonic",strain="MT72"),
                  data.frame(gene=b.31$X,l2fc=b.31$log2FoldChange,pvalue=b.31$padj,condition="Biofilm",strain="MT31"),
                  data.frame(gene=p.31$X,l2fc=p.31$log2FoldChange,pvalue=p.31$padj,condition="Planktonic",strain="MT31"),
                  data.frame(gene=b.55$X,l2fc=b.55$log2FoldChange,pvalue=b.55$padj,condition="Biofilm",strain="MT55"),
                  data.frame(gene=p.55$X,l2fc=p.55$log2FoldChange,pvalue=p.55$padj,condition="Planktonic",strain="MT55"))

all.lpda <- all.lpda[all.lpda$gene %in% operon.genes,]
all.lpda <- all.lpda %>% mutate(sig=case_when(pvalue < 0.05 ~ "Yes",
                                              pvalue >= 0.05 ~ "No" ))
all.lpda$sc <- paste(all.lpda$strain,all.lpda$condition,sep="-")
all.lpda$sc <- factor(all.lpda$sc,
             levels=c("MT31-Planktonic","MT31-Biofilm",
                      "MT345-Planktonic","MT345-Biofilm",
                      "MT49-Planktonic","MT49-Biofilm",
                      "MT55-Planktonic","MT55-Biofilm",
                      "MT72-Planktonic","MT72-Biofilm",
                      "MT540-Planktonic","MT540-Biofilm"
                      ))

strain.col <- c("#FFCCFF","#CC99FF","lightblue","dodgerblue","#FFCC99","#FF6666")
strain.order <- c("MT49","MT540","MT31","MT55","MT72","MT345")

lpda.l2fc <- ggplot(all.lpda[all.lpda$condition == "Planktonic",],aes(x=gene,y=l2fc)) + geom_point(aes(color=sig,fill=strain),size=4,pch=21,stroke=1)+
  facet_wrap(~sc)+ scale_x_discrete(limits=operon.genes)+
  theme_bw()+geom_hline(yintercept=0,linetype="dashed",size=0.25)+
  theme(axis.text.x=element_text(angle=30))+
  scale_fill_manual(name="Strain",values=strain.col,breaks=strain.order,labels=strain.order)+
  scale_color_manual(name="Significant?",values=c("black","white"),breaks=c("y","n"))+
  ylab("Log2 Fold Change")+xlab("Gene")

lpda.l2fc

ExportPlot(lpda.l2fc,"DE_lpdA/Figures/allStrains_lpdAoperon_planktonic",width=8,height=6)

lpda.l2fc <- ggplot(all.operon.counts[all.operon.counts$condition == "biofilm",],aes(x=gene,y=l2fc)) + geom_point(aes(color=sig,fill=strain),size=4,pch=21,stroke=1)+
  facet_wrap(~sc)+ scale_x_discrete(limits=operon.genes)+
  theme_bw()+geom_hline(yintercept=0,linetype="dashed",size=0.25)+
  theme(axis.text.x=element_text(angle=30))+
  scale_fill_manual(name="Strain",values=strain.col,breaks=strain.order,labels=strain.order)+
  scale_color_manual(name="Significant?",values=c("black","white"),breaks=c("y","n"))+
  ylab("Log2 Fold Change")+xlab("Gene")

lpda.l2fc

ExportPlot(lpda.l2fc,"DE_lpdA/Figures/allStrains_lpdAoperon_biofilm",width=8,height=6)

#strain.col <- c("#CC99FF","#FFCCFF")
strain.col <- c("lightblue","dodgerblue","#FF6666","#FFCC99","#FFCCFF","#CC99FF")

lpda.l2fc <- ggplot(all.lpda,aes(x=gene,y=l2fc)) + geom_point(aes(color=sig,fill=strain),size=4,pch=21,stroke=1)+
  facet_wrap(~sc)+ scale_x_discrete(limits=operon.genes)+
  theme_bw()+geom_hline(yintercept=0,linetype="dashed",size=0.25)+
  theme(axis.text.x=element_text(angle=30))+
  scale_fill_manual(name="strain",values=strain.col)+
  scale_color_manual(name="significant?",values=c("black","white"),breaks=c("Yes","No"))+
  ylab("log2 fold change") + xlab("Gene")+
  scale_y_continuous(limits=c(-2,4),breaks=c(-1,-2,0,1,2,3,4))
lpda.l2fc

ExportPlot(lpda.l2fc,"../NewFigures/Supplement/FigureS6",width=6,height=6)

# all strains
operon.genes <- c("lpdA","glpD2","phoY1","Rv3300c","atsB","lpqC")
b.540 <- read.csv("DEFiles/BEvA/540_BFEvA_allGenes.csv")
p.540 <- read.csv("DEFiles/PEvA/540_PEvA_allGenes.csv")
b.49 <- read.csv("DEFiles/BEvA/49_BFEvA_allGenes.csv")
p.49 <- read.csv("DEFiles/PEvA/49_PEvA_allGenes.csv")

b.345 <- read.csv("DEFiles/BEvA/345_BFEvA_allGenes.csv")
p.345 <- read.csv("DEFiles/PEvA/345_PEvA_allGenes.csv")
b.72 <- read.csv("DEFiles/BEvA/72_BFEvA_allGenes.csv")
p.72 <- read.csv("DEFiles/PEvA/72_PEvA_allGenes.csv")

b.31 <- read.csv("DEFiles/BEvA/31_BFEvA_allGenes.csv")
p.31 <- read.csv("DEFiles/PEvA/31_PEvA_allGenes.csv")
b.55 <- read.csv("DEFiles/BEvA/55_BFEvA_allGenes.csv")
p.55 <- read.csv("DEFiles/PEvA/55_PEvA_allGenes.csv")

all.lpda <- rbind(data.frame(gene=b.540$X,l2fc=b.540$log2FoldChange,pvalue=b.540$padj,condition="Biofilm",strain="MT540"),
                  data.frame(gene=p.540$X,l2fc=p.540$log2FoldChange,pvalue=p.540$padj,condition="Planktonic",strain="MT540"),
                  data.frame(gene=b.49$X,l2fc=b.49$log2FoldChange,pvalue=b.49$padj,condition="Biofilm",strain="MT49"),
                  data.frame(gene=p.49$X,l2fc=p.49$log2FoldChange,pvalue=p.49$padj,condition="Planktonic",strain="MT49"),
                  data.frame(gene=b.345$X,l2fc=b.345$log2FoldChange,pvalue=b.345$padj,condition="Biofilm",strain="MT345"),
                  data.frame(gene=p.345$X,l2fc=p.345$log2FoldChange,pvalue=p.345$padj,condition="Planktonic",strain="MT345"),
                  data.frame(gene=b.72$X,l2fc=b.72$log2FoldChange,pvalue=b.72$padj,condition="Biofilm",strain="MT72"),
                  data.frame(gene=p.72$X,l2fc=p.72$log2FoldChange,pvalue=p.72$padj,condition="Planktonic",strain="MT72"),
                  data.frame(gene=b.31$X,l2fc=b.31$log2FoldChange,pvalue=b.31$padj,condition="Biofilm",strain="MT31"),
                  data.frame(gene=p.31$X,l2fc=p.31$log2FoldChange,pvalue=p.31$padj,condition="Planktonic",strain="MT31"),
                  data.frame(gene=b.55$X,l2fc=b.55$log2FoldChange,pvalue=b.55$padj,condition="Biofilm",strain="MT55"),
                  data.frame(gene=p.55$X,l2fc=p.55$log2FoldChange,pvalue=p.55$padj,condition="Planktonic",strain="MT55"))

all.lpda <- all.lpda[all.lpda$gene %in% operon.genes,]
all.lpda <- all.lpda %>% mutate(sig=case_when(pvalue < 0.05 ~ "Yes",
                                              pvalue >= 0.05 ~ "No" ))
all.lpda$sc <- paste(all.lpda$strain,all.lpda$condition,sep="-")
all.lpda$sc <- factor(all.lpda$sc,
                      levels=c("MT31-Planktonic","MT31-Biofilm",
                               "MT345-Planktonic","MT345-Biofilm",
                               "MT49-Planktonic","MT49-Biofilm",
                               "MT55-Planktonic","MT55-Biofilm",
                               "MT72-Planktonic","MT72-Biofilm",
                               "MT540-Planktonic","MT540-Biofilm"
                      ))

strain.col <- c("#FFCCFF","#CC99FF","lightblue","dodgerblue","#FFCC99","#FF6666")
strain.order <- c("MT49","MT540","MT31","MT55","MT72","MT345")

lpda.l2fc <- ggplot(all.lpda[all.lpda$condition == "Planktonic",],aes(x=gene,y=l2fc)) + geom_point(aes(color=sig,fill=strain),size=4,pch=21,stroke=1)+
  facet_wrap(~sc)+ scale_x_discrete(limits=operon.genes)+
  theme_bw()+geom_hline(yintercept=0,linetype="dashed",size=0.25)+
  theme(axis.text.x=element_text(angle=30))+
  scale_fill_manual(name="Strain",values=strain.col,breaks=strain.order,labels=strain.order)+
  scale_color_manual(name="Significant?",values=c("black","white"),breaks=c("y","n"))+
  ylab("Log2 Fold Change")+xlab("Gene")

lpda.l2fc

ExportPlot(lpda.l2fc,"DE_lpdA/Figures/allStrains_lpdAoperon_planktonic",width=8,height=6)

lpda.l2fc <- ggplot(all.operon.counts[all.operon.counts$condition == "biofilm",],aes(x=gene,y=l2fc)) + geom_point(aes(color=sig,fill=strain),size=4,pch=21,stroke=1)+
  facet_wrap(~sc)+ scale_x_discrete(limits=operon.genes)+
  theme_bw()+geom_hline(yintercept=0,linetype="dashed",size=0.25)+
  theme(axis.text.x=element_text(angle=30))+
  scale_fill_manual(name="Strain",values=strain.col,breaks=strain.order,labels=strain.order)+
  scale_color_manual(name="Significant?",values=c("black","white"),breaks=c("y","n"))+
  ylab("Log2 Fold Change")+xlab("Gene")

lpda.l2fc

ExportPlot(lpda.l2fc,"DE_lpdA/Figures/allStrains_lpdAoperon_biofilm",width=8,height=6)
