require(DESeq2)
require(dplyr)
require(ggplot2)
require(RColorBrewer)
require(pheatmap)
require(viridis)
require(reshape2)
require(ggpubr)

#####
# Analysis of M. tb RNA-seq data from evolved populations, using DESeq2 and counts from HTSeq.
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
sampleData <- read.delim("Sautons_Stationary_metadata.txt", sep="", header = TRUE)

# reformat metadata
sampleTable <- data.frame(SampleName = sampleData$LibraryName,
                          CountsFile = sampleData$CountsFile,
                          Strain = factor(sampleData$Strain),
                          Genotype = factor(sampleData$Genotype),
                          Condition = factor(sampleData$Condition,levels = c("Planktonic","Biofilm")),
                          Phase = factor(sampleData$GrowthPhase, levels = c("Exponential","Stationary","Biofilm","Planktonic")),
                          Media = factor(sampleData$Media),
                          Clade = factor(sampleData$Clade),
                          SampleGroup = sampleData$SampleGroup,
                          SampleID = sampleData$SampleID)

### H37Rv analyses ###

# colors, etc
custom.palette <- colorRampPalette(rev(c("red","white","blue")))(20)
my_colors <- c("grey",custom.palette[15],custom.palette[5])

# all Rv PCA (Figure S8A)
t0 <- sampleTable[sampleTable$Strain == "H37Rv",]
t0.DESeq <- DESeqDataSetFromHTSeqCount(sampleTable = t0, directory = ".", design = ~ Phase)
t0.DESeq <- t0.DESeq[rowSums(counts(t0.DESeq)) > 1,]
t0.DESeq <- DESeq(t0.DESeq)
t0.vst <- vst(t0.DESeq, blind=FALSE)

t0.pca <- DESeq2::plotPCA(t0.vst, intgroup=c("Media","Phase")) + theme_bw()+
  theme(axis.text = element_text(size=10),axis.title = element_text(size=12), legend.text = element_text(size=10),
        legend.title = element_blank(),legend.position="right")
t0.pca

ExportPlot(t0.pca,"Figures/H37Rv_PCA",width=6,height=4)

# compare stationary and exponential phases in 7h9 planktonic culture
t1 <- sampleTable[sampleTable$Strain == "H37Rv" & sampleTable$Media == "7H9",]
t1.DESeq <- DESeqDataSetFromHTSeqCount(sampleTable = t1, directory = ".", design = ~ Phase)
t1.DESeq <- t1.DESeq[rowSums(counts(t1.DESeq)) > 1,]
t1.DESeq <- DESeq(t1.DESeq)
t1.result <- results(t1.DESeq, alpha = 0.05,lfcThreshold = 0, contrast=c("Phase","Stationary","Exponential"))
sum(t1.result$padj <0.05, na.rm=TRUE)
t1.sig <- data.frame(subset(t1.result, padj<0.05))
voldata.1 <- data.frame(t1.result)


# compare biofilm and planktonic conditions in sauton's media
t2 <- sampleTable[sampleTable$Strain == "H37Rv" & sampleTable$Media == "Sautons",]
t2.DESeq <- DESeqDataSetFromHTSeqCount(sampleTable = t2, directory = ".", design = ~ Condition)
t2.DESeq <- t2.DESeq[rowSums(counts(t2.DESeq)) > 1,]
t2.DESeq <- DESeq(t2.DESeq)
t2.result <- results(t2.DESeq, alpha = 0.05,lfcThreshold = 0, contrast=c("Condition","Biofilm","Planktonic"))
sum(t2.result$padj <0.05, na.rm=TRUE)
t2.sig <- data.frame(subset(t2.result, padj<0.05))
voldata.2 <- data.frame(t2.result)

# compare biofilm and planktonic conditions in sauton's and 7h9 media, respectively
t3 <- sampleTable[sampleTable$SampleName %in% c("H37Rv-HE-A","H37Rv-HE-B","H37Rv-HE-C","H37Rv-A","H37Rv-B","H37Rv-C"),]
t3.DESeq <- DESeqDataSetFromHTSeqCount(sampleTable = t3, directory = ".", design = ~ Condition)
t3.DESeq <- t3.DESeq[rowSums(counts(t3.DESeq)) > 1,]
t3.DESeq <- DESeq(t3.DESeq)
t3.result <- results(t3.DESeq, alpha = 0.05,lfcThreshold = 0, contrast=c("Condition","Biofilm","Planktonic"))
t3.sig <- data.frame(subset(t3.result, padj<0.05))
voldata.3 <- data.frame(t3.result)

# compare H37Rv planktonic culture in sauton's and 7h9 media
t4 <- sampleTable[sampleTable$Strain == "H37Rv" & sampleTable$Condition == "Planktonic",]
t4.DESeq <- DESeqDataSetFromHTSeqCount(sampleTable = t4, directory = ".", design = ~ Media)
t4.DESeq <- t4.DESeq[rowSums(counts(t4.DESeq)) > 1,]
t4.DESeq <- DESeq(t4.DESeq)
t4.result <- results(t4.DESeq, alpha = 0.05,lfcThreshold = 0, contrast=c("Media","Sautons","7H9"))
t4.sig <- data.frame(subset(t4.result, padj<0.05))
voldata.4 <- data.frame(t4.result)

### Pellicle evolved 31 analyses ###

# compare stationary and exponential phases in 7h9 planktonic culture
t5 <- sampleTable[sampleTable$Strain == "31" & sampleTable$Media == "7H9",]
t5.DESeq <- DESeqDataSetFromHTSeqCount(sampleTable = t5, directory = ".", design = ~ Phase)
t5.DESeq <- t5.DESeq[rowSums(counts(t5.DESeq)) > 1,]
t5.DESeq <- DESeq(t5.DESeq)
t5.result <- results(t5.DESeq, alpha = 0.05,lfcThreshold = 0, contrast=c("Phase","Stationary","Exponential"))
sum(t5.result$padj <0.05, na.rm=TRUE)
voldata.5 <- data.frame(t5.result)

# are stationary genes the same as biofilm genes? (evolved MT31)

dat.31 <- read.csv("../mtb_ExpEvo_RNA/DEFiles/EBvP/31_evoBvP.csv")
stat <- rbind(data.frame(type="stat",padj=voldata.5$padj,lfc=voldata.5$log2FoldChange,gene=rownames(voldata.5)),
              data.frame(type="bf",padj=dat.31$padj,lfc=dat.31$log2FoldChange,gene=dat.31$X))
stat.up <- stat$gene[stat$type == "stat" & stat$padj < 0.05 & stat$lfc > 2]; length(stat.up)
stat.down <- stat$gene[stat$type == "stat" & stat$padj < 0.05 & stat$lfc < -2]; length(stat.down)
bf.up <- stat$gene[stat$type == "bf" & stat$padj < 0.05 & stat$lfc > 2]; length(bf.up)
bf.down <- stat$gene[stat$type == "bf" & stat$padj < 0.05 & stat$lfc < -2]; length(bf.down)

up <- stat.up[stat.up %in% bf.up]; length(up)
down <- stat.down[stat.down %in% bf.down]; length(down)
stat.up.unique <- stat.up[!(stat.up %in% bf.up)]; length(stat.up.unique)
stat.down.unique <- stat.down[!(stat.down %in% bf.down)]; length(stat.down.unique)
bf.up.unique <- bf.up[!(bf.up %in% stat.up)]; length(bf.up.unique)
bf.down.unique <- bf.down[!(bf.down %in% stat.down)]; length(bf.down.unique)

mt31.bar.data <- data.frame(dir=rep(c("up","down"),3), lin=c("Stationary","Stationary","Biofilm","Biofilm","Biofilm & Stationary","Biofilm & Stationary"),
                       count=c(length(stat.up.unique),length(stat.down.unique),length(bf.up.unique),length(bf.down.unique),length(up),length(down)))
mt31.bar.data$lin <- factor(mt31.bar.data$lin,levels=c("Biofilm","Stationary","Biofilm & Stationary"))

mt31.count.plot <- ggplot(mt31.bar.data,(aes(x=lin,y=count))) + geom_col(aes(fill=lin,alpha=dir),width=0.5)+ theme_bw()+
  scale_alpha_manual(name=NULL,values=c(0.4,1),labels=c("downregulated","upregulated"),breaks=c("down","up"))+
  xlab(NULL)+ylab("Count")+
  scale_fill_manual(name=NULL,values=viridis(n=5,option="B",begin=0.2,end=0.8)[c(1,3,5)])+
  guides(fill=guide_legend(order=1))+
  theme(axis.text=element_text(size=12),axis.title.y = element_text(size=12),axis.text.x=element_text(angle=20))
mt31.count.plot

# biofilm vs sautons vs stationary in Rv
stat <- rbind(data.frame(type="Biofilm",padj=voldata.3$padj,lfc=voldata.3$log2FoldChange,gene=rownames(voldata.3)),
              data.frame(type="Sautons",padj=voldata.4$padj,lfc=voldata.4$log2FoldChange,gene=rownames(voldata.4)),
              data.frame(type="Stationary",padj=voldata.1$padj,lfc=voldata.1$log2FoldChange,gene=rownames(voldata.1)))

saut.up <- stat$gene[stat$type == "Sautons" & stat$padj < 0.05 & stat$lfc > 2]; length(saut.up)
saut.down <- stat$gene[stat$type == "Sautons" & stat$padj < 0.05 & stat$lfc < -2]; length(saut.down)
bf.up <- stat$gene[stat$type == "Biofilm" & stat$padj < 0.05 & stat$lfc > 2]; length(bf.up)
bf.down <- stat$gene[stat$type == "Biofilm" & stat$padj < 0.05 & stat$lfc < -2]; length(bf.down)
stat.up <- stat$gene[stat$type == "Stationary" & stat$padj < 0.05 & stat$lfc > 2]; length(stat.up)
stat.down <- stat$gene[stat$type == "Stationary" & stat$padj < 0.05 & stat$lfc < -2]; length(stat.down)

bf.saut.up <- saut.up[saut.up %in% bf.up]; length(bf.saut.up)
bf.saut.down <- saut.down[saut.down %in% bf.down]; length(bf.saut.down)
bf.stat.up <- stat.up[stat.up %in% bf.up]; length(bf.stat.up)
bf.stat.down <- stat.down[stat.down %in% bf.down]; length(bf.stat.down)
saut.up.unique <- saut.up[!(saut.up %in% bf.up)]; length(saut.up.unique)
saut.down.unique <- saut.down[!(saut.down %in% bf.down)]; length(saut.down.unique)
stat.up.unique <- stat.up[!(stat.up %in% bf.up)]; length(stat.up.unique)
stat.down.unique <- stat.down[!(stat.down %in% bf.down)]; length(stat.down.unique)
bf.up.unique <- bf.up[!(bf.up %in% stat.up)]; length(bf.up.unique)
bf.down.unique <- bf.down[!(bf.down %in% stat.down)]; length(bf.down.unique)

shared.genes <- rbind(data.frame(type="bf-saut",dir="up",gene=bf.saut.up),
                      data.frame(type="bf-saut",dir="down",gene=bf.saut.down),
                      data.frame(type="bf-stat",dir="up",gene=bf.stat.up),
                      data.frame(type="bf-stat",dir="down",gene=bf.stat.down))

rv.bar.data <- data.frame(dir=rep(c("up","down"),5),
                          lin=c(rep("Stationary",2),rep("Biofilm",2),rep("Sautons",2),
                                rep("Biofilm & Sautons",2),rep("Biofilm & Stationary",2)),
                          count=c(length(stat.up.unique),length(stat.down.unique),
                                  length(bf.up.unique),length(bf.down.unique),
                                  length(saut.up.unique),length(saut.down.unique),
                                  length(bf.saut.up),length(bf.saut.down),
                                  length(bf.stat.up),length(bf.stat.down)))

rv.bar.data$lin <- factor(rv.bar.data$lin, levels=c("Biofilm","Sautons","Stationary","Biofilm & Sautons","Biofilm & Stationary"))

rv.count.plot <- ggplot(rv.bar.data,(aes(x=lin,y=count))) + geom_col(aes(fill=lin,alpha=dir),width=0.5)+ theme_bw()+
  scale_alpha_manual(name=NULL,values=c(0.4,1),labels=c("downregulated","upregulated"),breaks=c("down","up"))+
  xlab(NULL)+ylab("Count")+
  scale_fill_manual(name=NULL,values=viridis(n=5,option="B",begin=0.2,end=0.8))+
  guides(fill=guide_legend(order=1))+
  theme(axis.text=element_text(size=12),axis.title.y = element_text(size=12),axis.text.x = element_text(angle=20))
rv.count.plot


all.bar.data <- rbind(data.frame(lin=rv.bar.data$lin,dir=rv.bar.data$dir,count=rv.bar.data$count,strain="H37Rv"),
                      data.frame(lin=mt31.bar.data$lin,dir=mt31.bar.data$dir,count=mt31.bar.data$count,strain="MT31-12"))

# Figure S8B
all.count.plot <- ggplot(all.bar.data,(aes(x=lin,y=count))) + geom_col(aes(fill=lin,alpha=dir),width=0.5)+ theme_bw()+
  scale_alpha_manual(name=NULL,values=c(0.4,1),labels=c("downregulated","upregulated"),breaks=c("down","up"))+
  xlab(NULL)+ylab("Count")+
  scale_fill_manual(name=NULL,values=viridis(n=5,option="B",begin=0.2,end=0.8))+
  guides(fill=guide_legend(order=1))+
  theme(axis.text=element_text(size=12),axis.title.y = element_text(size=12),axis.text.x = element_text(angle=20))+
  facet_wrap(~strain, drop=T, scales="free_x")
all.count.plot
