require(DESeq2)
require(dplyr)
require(ggplot2)
require(RColorBrewer)
require(pheatmap)
require(viridis)
require(reshape2)
require(ggpubr)
require(scales)

#####
# Differential expression analysis of duplicated genes, using DESeq2 and counts from HTSeq.
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

#### SK31 ####
dupgenes.31 <- read.delim("Metadata/31_DupGenes.txt",header=F)$V1
B31.table <- sampleTable[sampleTable$Strain == "31" & sampleTable$Condition == "Biofilm",]

B31.DESeq <- DESeqDataSetFromHTSeqCount(sampleTable = B31.table, directory = ".", design = ~ Genotype)
B31.DESeq <- B31.DESeq[rowSums(counts(B31.DESeq)) > 1,]
B31.DESeq <- DESeq(B31.DESeq)

r.31 <- results(B31.DESeq, alpha = 0.05, lfcThreshold = 0, contrast=c("Genotype","Evolved","Ancestral"))

r.31.dup <- subset(r.31, rownames(r.31) %in% dupgenes.31)
r.31.non.dup <- subset(r.31, rownames(r.31) %in% setdiff(rownames(assay(B31.DESeq)),dupgenes.31))

#write.csv(x = r.31.dup, file="DEFiles/Duplication/31_dupGenes_l2fc.csv")
#write.csv(x= r.31.non.dup, file="DEFiles/Duplication/31_nonDup_l2fc.csv")

#### SK55 ####
dupgenes.55 <- read.delim("Metadata/55_DupGenes.txt",header=F)$V1

B55.table <- sampleTable[sampleTable$Strain == "55" & sampleTable$Condition == "Biofilm",]

B55.DESeq <- DESeqDataSetFromHTSeqCount(sampleTable = B55.table, directory = ".", design = ~ Genotype)
B55.DESeq <- B55.DESeq[rowSums(counts(B55.DESeq)) > 1,]
B55.DESeq <- DESeq(B55.DESeq)

r.55 <- results(B55.DESeq, alpha = 0.05, lfcThreshold = 0, contrast=c("Genotype","Evolved","Ancestral"))

r.55.dup <- subset(r.55, rownames(r.55) %in% dupgenes.55)
r.55.non.dup <- subset(r.55, rownames(r.55) %in% setdiff(rownames(assay(B55.DESeq)),dupgenes.55))

#write.csv(x = r.55.dup, file="DEFiles/Duplication/55_dupGenes_l2fc.csv")
#write.csv(x= r.55.non.dup, file="DEFiles/Duplication/55_nonDup_l2fc.csv")

#### plot l2fc ####

dupcounts.31 <- read.csv("DEFiles/Duplication/31_dupGenes_l2fc.csv",header=T)
dupcounts.55 <- read.csv("DEFiles/Duplication/55_dupGenes_l2fc.csv",header=T)
counts.31 <- read.csv("DEFiles/Duplication/31_nonDup_l2fc.csv",header=T)
counts.55 <- read.csv("DEFiles/Duplication/55_nonDup_l2fc.csv", header=T)

counts31 <- data.frame(count=counts.31$log2FoldChange,gene=counts.31$X,strain="31")
dupcounts31 <- data.frame(count=dupcounts.31$log2FoldChange,gene=dupcounts.31$X,strain="31")
counts55 <- data.frame(count=counts.55$log2FoldChange,gene=counts.55$X,strain="55")
dupcounts55 <- data.frame(count=dupcounts.55$log2FoldChange,gene=dupcounts.55$X,strain="55")

allcounts <- rbind(data.frame(gene=counts31$gene,l2fc=counts31$count,type="all",strain="31"),
                   data.frame(gene=dupcounts31$gene,l2fc=dupcounts31$count,type="duplication",strain="31"),
                   data.frame(gene=counts55$gene,l2fc=counts55$count,type="all",strain="55"),
                   data.frame(gene=dupcounts55$gene,l2fc=dupcounts55$count,type="duplication",strain="55"))

allcounts$strain <- as.factor(allcounts$strain)

# boxplots comparing l2fc values of duplicated vs non-duplicated genes (Figure 6C-right)
box31 <- ggplot(allcounts[allcounts$strain=="31",],aes(x=type,y=l2fc)) + geom_hline(yintercept=0,linetype="dashed")+
  geom_boxplot(alpha=1,aes(fill=type)) + theme_minimal()+
  theme(legend.position="none",axis.title = element_blank(),
        axis.text.y=element_text(size=10),axis.text.x=element_blank())+
  scale_fill_manual(values=c("lightblue","darkblue"),breaks=c("all","duplication"))+
  scale_y_continuous(limits=c(-5,5),breaks=c(-5,-2.5,0,2.5,5))
  
box31

stat31 <- compare_means(l2fc ~ type, data = allcounts[allcounts$strain=="31",], p.adjust.method = "BH")
stat31 <- stat31 %>% mutate(y.position=c(1.5))
box31.stats <- box31 + stat_pvalue_manual(stat31, label="p.signif",label.size=5)
box31.stats

box55 <- ggplot(allcounts[allcounts$strain=="55",],aes(x=type,y=l2fc)) +geom_hline(yintercept=0,linetype="dashed")+
  geom_boxplot(alpha=1,aes(fill=type)) + theme_minimal()+
  theme(legend.position="none",axis.title = element_blank(),
        axis.text.y=element_text(size=10),axis.text.x=element_blank())+
  scale_fill_manual(values=c("lightblue","darkblue"),breaks=c("all","duplication"))+
  scale_y_continuous(limits=c(-10,10),breaks=c(-10,-5,0,5,10))
box55

stat55 <- compare_means(l2fc ~ type, data = allcounts[allcounts$strain=="55",], p.adjust.method = "BH")
stat55 <- stat55 %>% mutate(y.position=c(4))
box55.stats <- box55 + stat_pvalue_manual(stat55, label="p.signif",label.size=5)
box55.stats

# sliding window plots (Figure 6C-left)
genes.pos <- read.delim("Metadata/genepos.txt",header=F)
data31 <- allcounts[allcounts$strain == "31",]

genes.pos <- genes.pos[-c(which(genes.pos$V1 %in% setdiff(genes.pos$V1,data31$gene))),]
data31$pos[match(genes.pos$V1,data31$gene)] <- genes.pos$V2

plot31 <- ggplot(data31,aes(x=pos,y=l2fc)) + geom_point(aes(color=type))+
  theme(axis.text = element_text(size=16),
        axis.title=element_text(size=14),
        legend.text = element_text(size=16),legend.title=element_blank())+
  scale_color_manual(name=NULL,values=c("lightblue","navy"),breaks=c("all","duplication"),
                     labels=c("Not duplicated","Duplicated"))+
  geom_hline(yintercept=0,linetype="dashed")+
  xlab(NULL)+ylab("log2 fold change")+theme_minimal()+
  scale_x_continuous(labels=comma)+
  scale_y_continuous(limits=c(-5,5),breaks=c(-5,-2.5,0,2.5,5))

plot31


data55 <- allcounts[allcounts$strain == "55",]
genes.pos <- genes.pos[-c(which(genes.pos$V1 %in% setdiff(genes.pos$V1,data55$gene))),]
data55$pos[match(genes.pos$V1,data55$gene)] <- genes.pos$V2


plot55 <- ggplot(data55,aes(x=pos,y=l2fc)) + geom_point(aes(color=type))+
  theme(axis.text = element_text(size=16),
        axis.title=element_text(size=14),
        legend.text = element_text(size=16),legend.title=element_blank())+
  scale_color_manual(name=NULL,values=c("lightblue","navy"),breaks=c("all","duplication"),
                     labels=c("Not duplicated","Duplicated"))+
  geom_hline(yintercept=0,linetype="dashed")+
  xlab("Genome position")+ylab("log2 fold change")+
  theme_minimal()+scale_x_continuous(labels=comma)+
  scale_y_continuous(limits=c(-10,10),breaks=c(-10,-5,0,5,10))

plot55


full.fig <- ggarrange(plot31,box31.stats,plot55,box55.stats,labels=c("MT31",NA,"MT55",NA),nrow=2,ncol=2,common.legend = T,
                      widths=c(3,1))
full.fig

# comparing duplicated to non-duplicated genes for other strains
strains <- c("31","55","49","540","72","345")
dupgenes <- read.delim("Metadata/dupgenes.txt",header=F)$V1

for (strain in strains){
  bst <- sampleTable[sampleTable$Strain == strain & sampleTable$Condition == "Biofilm",]
  pst <- sampleTable[sampleTable$Strain == strain & sampleTable$Condition == "Planktonic",]
  bt <- DESeqDataSetFromHTSeqCount(sampleTable = bst, directory = ".", design = ~ Genotype)
  bt <- bt[rowSums(counts(bt)) > 1,]
  bt <- DESeq(bt)
  pt <- DESeqDataSetFromHTSeqCount(sampleTable = pst, directory = ".", design = ~ Genotype)
  pt <- pt[rowSums(counts(pt)) > 1,]
  pt <- DESeq(pt)
  br <- results(bt, alpha = 0.05, contrast=c("Genotype","Evolved","Ancestral"))
  pr <- results(pt, alpha = 0.05, contrast=c("Genotype","Evolved","Ancestral"))
  bdup <- subset(br,rownames(br) %in% dupgenes)
  bnondup <- subset(br, rownames(br) %in% setdiff(rownames(assay(bt)),dupgenes))
  pdup <- subset(pr,rownames(pr) %in% dupgenes)
  pnondup <- subset(pr, rownames(pr) %in% setdiff(rownames(assay(pt)),dupgenes))
  write.csv(x=bdup,file=paste("DEFiles/Duplication/",strain,"_BdupGenes.csv",sep=""))
  write.csv(x=bnondup,file=paste("DEFiles/Duplication/",strain,"_BnonDup.csv",sep=""))
  write.csv(x=pdup,file=paste("DEFiles/Duplication/",strain,"_PdupGenes.csv",sep=""))
  write.csv(x=pnondup,file=paste("DEFiles/Duplication/",strain,"_PnonDup.csv",sep=""))
}

# plot l2fc values from biofilm samples (Figure S9A)

dupcounts.31 <- read.csv("DEFiles/Duplication/31_BdupGenes.csv",header=T)
counts.31 <- read.csv("DEFiles/Duplication/31_BnonDup.csv",header=T)
dupcounts.55 <- read.csv("DEFiles/Duplication/55_BdupGenes.csv",header=T)
counts.55 <- read.csv("DEFiles/Duplication/55_BnonDup.csv", header=T)
dupcounts.49 <- read.csv("DEFiles/Duplication/49_BdupGenes.csv",header=T)
counts.49 <- read.csv("DEFiles/Duplication/49_BnonDup.csv",header=T)
dupcounts.540 <- read.csv("DEFiles/Duplication/540_BdupGenes.csv",header=T)
counts.540 <- read.csv("DEFiles/Duplication/540_BnonDup.csv", header=T)
dupcounts.72 <- read.csv("DEFiles/Duplication/72_BdupGenes.csv",header=T)
counts.72 <- read.csv("DEFiles/Duplication/72_BnonDup.csv",header=T)
dupcounts.345 <- read.csv("DEFiles/Duplication/345_BdupGenes.csv",header=T)
counts.345 <- read.csv("DEFiles/Duplication/345_BnonDup.csv", header=T)

allcounts <- rbind(data.frame(gene=counts.31$X,l2fc=counts.31$log2FoldChange,type="all",strain="MT31"),
                   data.frame(gene=dupcounts.31$X,l2fc=dupcounts.31$log2FoldChange,type="duplication",strain="MT31"),
                   data.frame(gene=counts.55$X,l2fc=counts.55$log2FoldChange,type="all",strain="MT55"),
                   data.frame(gene=dupcounts.55$X,l2fc=dupcounts.55$log2FoldChange,type="duplication",strain="MT55"),
                   data.frame(gene=counts.49$X,l2fc=counts.49$log2FoldChange,type="all",strain="MT49"),
                   data.frame(gene=dupcounts.49$X,l2fc=dupcounts.49$log2FoldChange,type="duplication",strain="MT49"),
                   data.frame(gene=counts.540$X,l2fc=counts.540$log2FoldChange,type="all",strain="MT540"),
                   data.frame(gene=dupcounts.540$X,l2fc=dupcounts.540$log2FoldChange,type="duplication",strain="MT540"),
                   data.frame(gene=counts.72$X,l2fc=counts.72$log2FoldChange,type="all",strain="MT72"),
                   data.frame(gene=dupcounts.72$X,l2fc=dupcounts.72$log2FoldChange,type="duplication",strain="MT72"),
                   data.frame(gene=counts.345$X,l2fc=counts.345$log2FoldChange,type="all",strain="MT345"),
                   data.frame(gene=dupcounts.345$X,l2fc=dupcounts.345$log2FoldChange,type="duplication",strain="MT345"))

strain.order <- c("MT31","MT49","MT72","MT55","MT540","MT345")
allcounts$strain <- factor(allcounts$strain,levels=strain.order)

all.box <- ggplot(allcounts,aes(x=type,y=l2fc)) + geom_hline(yintercept=0,linetype="dashed")+
  geom_boxplot(alpha=1,aes(fill=type)) + theme_bw()+
  facet_wrap(~strain,scales="free_y")+
  theme(legend.position="top",axis.title = element_blank(),
        axis.text.y=element_text(size=12),axis.text.x=element_blank())+
  scale_fill_manual(values=c("lightblue","darkblue"),breaks=c("all","duplication"))
all.box

ExportPlot(all.box,"../NewFigures/Supplement/FigureS4A",width=8,height=6)

# making individual boxplots to get the p-values
for (strain in strain.order){
  box <- ggplot(allcounts[allcounts$strain == strain,],aes(x=type,y=l2fc)) + geom_hline(yintercept=0,linetype="dashed")+
    geom_boxplot(alpha=1,aes(fill=type)) + theme_minimal()+
    theme(legend.position="none",axis.title = element_blank(),
          axis.text.y=element_text(size=10),axis.text.x=element_blank())+
    scale_fill_manual(values=c("lightblue","darkblue"),breaks=c("all","duplication"))
  stat <- compare_means(l2fc ~ type, data = allcounts[allcounts$strain==strain,], p.adjust.method = "BH")
  print(strain)
  print(stat)
  stat <- stat %>% mutate(y.position=c(5))
  assign(paste("box",strain,sep="."), box + stat_pvalue_manual(stat, label="p.adj",label.size=5))
}

# plot l2fc values from planktonic samples (Figure S9B)

dupcounts.31 <- read.csv("DEFiles/Duplication/31_PdupGenes.csv",header=T)
counts.31 <- read.csv("DEFiles/Duplication/31_PnonDup.csv",header=T)
dupcounts.55 <- read.csv("DEFiles/Duplication/55_PdupGenes.csv",header=T)
counts.55 <- read.csv("DEFiles/Duplication/55_PnonDup.csv", header=T)
dupcounts.49 <- read.csv("DEFiles/Duplication/49_PdupGenes.csv",header=T)
counts.49 <- read.csv("DEFiles/Duplication/49_PnonDup.csv",header=T)
dupcounts.540 <- read.csv("DEFiles/Duplication/540_PdupGenes.csv",header=T)
counts.540 <- read.csv("DEFiles/Duplication/540_PnonDup.csv", header=T)
dupcounts.72 <- read.csv("DEFiles/Duplication/72_PdupGenes.csv",header=T)
counts.72 <- read.csv("DEFiles/Duplication/72_PnonDup.csv",header=T)
dupcounts.345 <- read.csv("DEFiles/Duplication/345_PdupGenes.csv",header=T)
counts.345 <- read.csv("DEFiles/Duplication/345_PnonDup.csv", header=T)

allcounts <- rbind(data.frame(gene=counts.31$X,l2fc=counts.31$log2FoldChange,type="all",strain="MT31"),
                   data.frame(gene=dupcounts.31$X,l2fc=dupcounts.31$log2FoldChange,type="duplication",strain="MT31"),
                   data.frame(gene=counts.55$X,l2fc=counts.55$log2FoldChange,type="all",strain="MT55"),
                   data.frame(gene=dupcounts.55$X,l2fc=dupcounts.55$log2FoldChange,type="duplication",strain="MT55"),
                   data.frame(gene=counts.49$X,l2fc=counts.49$log2FoldChange,type="all",strain="MT49"),
                   data.frame(gene=dupcounts.49$X,l2fc=dupcounts.49$log2FoldChange,type="duplication",strain="MT49"),
                   data.frame(gene=counts.540$X,l2fc=counts.540$log2FoldChange,type="all",strain="MT540"),
                   data.frame(gene=dupcounts.540$X,l2fc=dupcounts.540$log2FoldChange,type="duplication",strain="MT540"),
                   data.frame(gene=counts.72$X,l2fc=counts.72$log2FoldChange,type="all",strain="MT72"),
                   data.frame(gene=dupcounts.72$X,l2fc=dupcounts.72$log2FoldChange,type="duplication",strain="MT72"),
                   data.frame(gene=counts.345$X,l2fc=counts.345$log2FoldChange,type="all",strain="MT345"),
                   data.frame(gene=dupcounts.345$X,l2fc=dupcounts.345$log2FoldChange,type="duplication",strain="MT345"))

strain.order <- c("MT31","MT49","MT72","MT55","MT540","MT345")
allcounts$strain <- factor(allcounts$strain,levels=strain.order)

all.box <- ggplot(allcounts,aes(x=type,y=l2fc)) + geom_hline(yintercept=0,linetype="dashed")+
  geom_boxplot(alpha=1,aes(fill=type)) + theme_bw()+
  facet_wrap(~strain,scales="free_y")+
  theme(legend.position="top",axis.title = element_blank(),
        axis.text.y=element_text(size=12),axis.text.x=element_blank())+
  scale_fill_manual(values=c("lightblue","darkblue"),breaks=c("all","duplication"))
all.box

# making individual boxplots to get the p-values
for (strain in strain.order){
  box <- ggplot(allcounts[allcounts$strain == strain,],aes(x=type,y=l2fc)) + geom_hline(yintercept=0,linetype="dashed")+
    geom_boxplot(alpha=1,aes(fill=type)) + theme_minimal()+
    theme(legend.position="none",axis.title = element_blank(),
          axis.text.y=element_text(size=10),axis.text.x=element_blank())+
    scale_fill_manual(values=c("lightblue","darkblue"),breaks=c("all","duplication"))
  stat <- compare_means(l2fc ~ type, data = allcounts[allcounts$strain==strain,], p.adjust.method = "BH")
  print(strain)
  print(stat)
  stat <- stat %>% mutate(y.position=c(5))
  assign(paste("box",strain,sep="."), box + stat_pvalue_manual(stat, label="p.adj",label.size=5))
}

# MT49 sliding window coverage (Figure S10A)
genes.pos <- read.delim("Metadata/genepos.txt",header=F)

dupcounts.49 <- read.csv("DEFiles/Duplication/49_BdupGenes.csv",header=T)
counts.49 <- read.csv("DEFiles/Duplication/49_BnonDup.csv",header=T)

data49 <- rbind(data.frame(gene=counts.49$X,l2fc=counts.49$log2FoldChange,type="all"),
                   data.frame(gene=dupcounts.49$X,l2fc=dupcounts.49$log2FoldChange,type="duplication"))

genes.pos <- genes.pos[-c(which(genes.pos$V1 %in% setdiff(genes.pos$V1,data49$gene))),]
data49$pos[match(genes.pos$V1,data49$gene)] <- genes.pos$V2

plot49 <- ggplot(data49,aes(x=pos,y=l2fc)) + geom_point(aes(color=type))+
  theme(axis.text.y = element_text(size=12),
        axis.title=element_text(size=14),
        legend.text = element_text(size=14),legend.title=element_blank())+
  scale_color_manual(name=NULL,values=c("lightblue","navy"),breaks=c("all","duplication"),
                     labels=c("Not duplicated","Duplicated"))+
  geom_hline(yintercept=0,linetype="dashed")+
  xlab("Genome position")+ylab("log2 fold change")+theme_minimal()+
  scale_x_continuous(labels=comma)

plot49

plot49.zoom <- plot49 + ylim(-2.5,2.5) + scale_x_continuous(limits=c(3250000,4000000),breaks = c(3250000,3500000,3750000,4000000)) +theme(legend.position="top")
plot49.zoom

# sliding window coverage for MT49 DNA and RNA seq (Figure S10B)
r.data <- read.delim("Metadata/49_RNAseq_coverage.txt",header=T)
d.data <- read.delim("Metadata/49_DNAseq_coverage.txt")[,c(1,2,4)]

all.data <- rbind(data.frame(position=r.data$position,coverage=r.data$X13.A.B,sample="X13.A.B",type="RNA",genotype="ancestral"),
                  data.frame(position=r.data$position,coverage=r.data$X14.A.B,sample="X14.A.B",type="RNA",genotype="ancestral"),
                  data.frame(position=r.data$position,coverage=r.data$X15.E.B,sample="X15.E.B",type="RNA",genotype="evolved"),
                  data.frame(position=r.data$position,coverage=r.data$X16.E.B,sample="X16.E.B",type="RNA",genotype="evolved"),
                  data.frame(position=d.data$position,coverage=d.data$X49.00,sample="X49.00",type="DNA",genotype="ancestral"),
                  data.frame(position=d.data$position,coverage=d.data$X49.12,sample="X49.12",type="DNA",genotype="evolved"))

plot <- ggplot(data = all.data, aes(x=position,y=coverage,color=sample)) + geom_line(linewidth=1)+
  ylab("Relative Sequencing Coverage\n")+theme_minimal()+xlab("\nGenome Position")+
  scale_color_manual(name="Genotype",breaks=c("X13.A.B","X14.A.B","X15.E.B","X16.E.B","X49.00","X49.12"),
                     values=c("#440154FF","#440154FF","#7AD151FF","#7AD151FF","#440154FF","#7AD151FF"))+
  scale_x_continuous(limits=c(3250000,4000000),breaks = c(3250000,3500000,3750000,4000000))+
  scale_y_continuous(limits=c(0,1.75),breaks=c(0,0.5,1.0,1.5))+
  theme(legend.position="top") + facet_wrap(~type,nrow=2)
plot

