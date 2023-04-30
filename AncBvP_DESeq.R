require(DESeq2)
require(dplyr)
require(ggplot2)
require(RColorBrewer)
require(pheatmap)
require(viridis)
require(reshape2)
require(tidyr)

#####
# Analysis of M. tb RNA-seq data from ancestral populations, using DESeq2 and counts from HTSeq.
# Requires: metadata file, HTSeq count files 
#####

setwd("~/Desktop/2023.04.26_RNAseq/mtb_ExpEvo_RNA/")

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

# subset to only ancestral samples
A.table <- sampleTable[sampleTable$Genotype == "Ancestral",]

# DE calculations by condition (biofilm vs plantkonic)
A.DESeq <- DESeqDataSetFromHTSeqCount(sampleTable = A.table, directory = ".", design = ~ Condition)

# filter out genes with 0 counts in all samples
A.DESeq <- A.DESeq[rowSums(counts(A.DESeq)) > 1,]

# run DESeq analysis
A.DESeq <- DESeq(A.DESeq)
resultsNames(A.DESeq)

# r-log transformation
A.rld <- rlogTransformation(A.DESeq, blind=FALSE)

# pca plot colored by growth condition (Figure 3B)
anc.pca.plot <- DESeq2::plotPCA(A.rld, intgroup=c("Condition")) + 
  theme_minimal()+ scale_color_manual(name="Growth Condition",values=c("#150E37FF","#D1426FFF"),
                                      breaks=c("Biofilm","Planktonic")) +
  theme(axis.text = element_text(size=10),axis.title = element_text(size=12), legend.text = element_text(size=12))
anc.pca.plot

# pca plot colored by strain (Figure S1A)
anc.pca.data <- DESeq2::plotPCA(A.rld, intgroup=c("Condition","Strain","WetWeight"),returnData=T) 
anc.pca.plot.2 <- ggplot(anc.pca.data, aes(x=PC1,y=PC2,color=Strain,shape=Condition)) + geom_point(size=4) +
  theme_minimal()+ scale_shape_manual(name="Growth Condition",values=c(16,17),breaks=c("Biofilm","Planktonic"))+
  xlab("PC1: 71% variance")+ylab("PC2: 13% variance") + scale_color_discrete(name="Population")+
  theme(axis.text = element_text(size=10),axis.title = element_text(size=12), 
        legend.text = element_text(size=12))
anc.pca.plot.2

ExportPlot(anc.pca.plot,"../NewFigures/Figure_3B",width=6,height=4)
ExportPlot(anc.pca.plot.2,"../NewFigures/Supplement/Figure_S1A",width=6,height=4)

# pca plot colored by wet weight (Figure S1B)
ww.pca.plot <-  ggplot(anc.pca.data,aes(x=PC1,y=PC2))+ theme_minimal()+
  geom_point(aes(color=WetWeight,shape=Condition),size=4)+
  scale_color_viridis(option = "A",end=0.9,begin=0.1,direction = -1)+
  scale_shape_manual(values=c(16,17),breaks=c("Biofilm","Planktonic"))+
  xlab("PC1: 71% variance")+ylab("PC2: 13% variance") +
  theme(axis.text = element_text(size=10),axis.title = element_text(size=12), 
        legend.text = element_text(size=12))
ww.pca.plot

ExportPlot(ww.pca.plot,"../NewFigures/Supplement/Figure_S1B",width=6,height=4)

# results of biofilm vs planktonic comparison
ABvP.result <- results(A.DESeq, alpha = 0.05,lfcThreshold = 0, contrast=c("Condition","Biofilm","Planktonic"))
write.csv(x = ABvP.result,"DEFiles/ABvP/all_ABvP.csv",quote=F,row.names = F)
# number of significantly DE genes
sum(ABvP.result$padj <0.05, na.rm=TRUE)
# subset all genes for only significant results (Supplementary Data 1 - AllPopulations_ABvP)
subset <- subset(ABvP.result,padj<0.05)
write.csv(x = subset,"DEFiles/ABvP/all_ABvP_DE.csv",quote=F,row.names = F)

# heatmap of all genes (Figure S2A)
my_colors <- magma(n=4, end = 0.8, begin=0.1)
condition <- my_colors[c(1,3)]
genotype <- my_colors[c(2,4)]
names(condition) <- c("Biofilm","Planktonic")
anno_colors <- list(Condition=condition)
BvPCounts <- assay(A.rld)[rownames(ABvP.result),]
BvPCounts <- BvPCounts - rowMeans(BvPCounts)
anno <- as.data.frame(colData(A.rld)[c("Condition")])
custom.palette <- c("red","white","blue")
colors <- colorRampPalette(rev(custom.palette))(20)

de.heatmap <- pheatmap(BvPCounts, color = colors, show_rownames = F, annotation_col = anno,show_colnames = F, 
                       cutree_cols = 1,treeheight_row = 0, annotation_colors = anno_colors,
                       kmeans_k = NA, border_color = "white", scale= "column", fontsize = 12)

ExportPlot(de.heatmap,"../NewFigures/Supplement/Figure_S2A",width=8,height=10)

#volcano plot (Figure 3C)
voldata <- data.frame(ABvP.result)
voldata <- voldata %>% mutate(sig=case_when(padj < 0.05 & log2FoldChange > 0 ~ "Significant-Up",
                                            padj < 0.05 & log2FoldChange < 0 ~ "Significant-Down",
                                            padj > 0.05 ~ "Not significant"))

custom.palette <- colorRampPalette(rev(c("red","white","blue")))(20)
my_colors <- c("grey",custom.palette[15],custom.palette[5])

vol <- ggplot(voldata,aes(x=log2FoldChange,y=-log10(padj))) + geom_point(aes(color=sig),alpha=0.8) + 
  theme_minimal()+xlim(-11,11)+ylim(0,80)+
  scale_color_manual(name=NULL,values=my_colors,breaks=c("Not significant","Significant-Up","Significant-Down"),
                     labels=c("Not significant","Upregulated","Downregulated"))+
  theme(legend.position = "top",axis.text = element_text(size=10),axis.title = element_text(size=12), 
        legend.text = element_text(size=12)) + xlab("log2 fold change")
vol

ExportPlot(vol,"../NewFigures/Figure_3C",width=6,height=6)


# separate DE analysis by individual strains (Supplementary Data 1 - ABvP for each strain)
strains <- c("31","49","55","72","345","540")

for (strain in strains){
  st <- A.table[A.table$Strain == strain,]
  t <- DESeqDataSetFromHTSeqCount(sampleTable = st, directory = ".", design = ~ Condition)
  t <- t[rowSums(counts(t)) > 1,]
  t <- DESeq(t)
  r <- results(t, alpha = 0.05, lfcThreshold = 0, contrast=c("Condition","Biofilm","Planktonic"))
  sub <- subset(r, padj < 0.05)
  write.csv(sub,paste("DEFiles/ABvP/",paste(strain,"ancBvP",sep="_"),".csv",sep=""),quote=F)
}

# reading individual data files
x31 <- data.frame(read.csv("DEFiles/ABvP/31_ancBvP.csv",header=T))
x49 <- data.frame(read.csv("DEFiles/ABvP/49_ancBvP.csv",header=T))
x55 <- data.frame(read.csv("DEFiles/ABvP/55_ancBvP.csv",header=T))
x72 <- data.frame(read.csv("DEFiles/ABvP/72_ancBvP.csv",header=T))
x345 <- data.frame(read.csv("DEFiles/ABvP/345_ancBvP.csv",header=T))
x540 <- data.frame(read.csv("DEFiles/ABvP/540_ancBvP.csv",header=T))

# reformatting data into one df
data <- rbind(data.frame(strain="31", gene = x31[x31$log2FoldChange >0,]$X, direction="UP"),
              data.frame(strain="31", gene = x31[x31$log2FoldChange <0,]$X, direction="DOWN"),
              data.frame(strain="49", gene = x49[x49$log2FoldChange >0,]$X, direction="UP"),
              data.frame(strain="49", gene = x49[x49$log2FoldChange <0,]$X, direction="DOWN"),
              data.frame(strain="55", gene = x55[x55$log2FoldChange >0,]$X, direction="UP"),
              data.frame(strain="55", gene = x55[x55$log2FoldChange <0,]$X, direction="DOWN"),
              data.frame(strain="72", gene = x72[x72$log2FoldChange >0,]$X, direction="UP"),
              data.frame(strain="72", gene = x72[x72$log2FoldChange <0,]$X, direction="DOWN"),
              data.frame(strain="345", gene = x345[x345$log2FoldChange >0,]$X, direction="UP"),
              data.frame(strain="345", gene = x345[x345$log2FoldChange <0,]$X, direction="DOWN"),
              data.frame(strain="540", gene = x540[x540$log2FoldChange >0,]$X, direction="UP"),
              data.frame(strain="540", gene = x540[x540$log2FoldChange <0,]$X, direction="DOWN"))

# color palette for all plots
custom.palette <- colorRampPalette(rev(c("red","white","blue")))(20)
my_colors <- c(custom.palette[5],custom.palette[15])

# correct order of strains for plotting
strains <- rev(c("31","55","345","72","49","540"))
strain.names <- rev(c("MT31","MT55","MT345","MT72","MT49","MT540"))
strain.order <- rev(c("x31","x55","x345","x72","x49","x540"))

# matrix of yes/no DE genes for each comparison
genes <- data.frame("gene"=unique(data$gene))

genes.2 <- genes %>% mutate("x31"=case_when(genes$gene %in% data[data$strain == "31" & data$direction == "UP",]$gene ~ 2,
                                            genes$gene %in% data[data$strain == "31" & data$direction == "DOWN",]$gene ~ 1),
                            "x49"=case_when(genes$gene %in% data[data$strain == "49" & data$direction == "UP",]$gene ~ 2,
                                            genes$gene %in% data[data$strain == "49" & data$direction == "DOWN",]$gene ~ 1),
                            "x55"=case_when(genes$gene %in% data[data$strain == "55" & data$direction == "UP",]$gene ~ 2,
                                            genes$gene %in% data[data$strain == "55" & data$direction == "DOWN",]$gene ~ 1),
                            "x72"=case_when(genes$gene %in% data[data$strain == "72" & data$direction == "UP",]$gene ~ 2,
                                            genes$gene %in% data[data$strain == "72" & data$direction == "DOWN",]$gene ~ 1),
                            "x345"=case_when(genes$gene %in% data[data$strain == "345" & data$direction == "UP",]$gene ~ 2,
                                             genes$gene %in% data[data$strain == "345" & data$direction == "DOWN",]$gene ~ 1),
                            "x540"=case_when(genes$gene %in% data[data$strain == "540" & data$direction == "UP",]$gene ~ 2,
                                             genes$gene %in% data[data$strain == "540" & data$direction == "DOWN",]$gene ~ 1))

genes.melt <- melt(genes.2,id.vars=c("gene"),measure.vars=strain.order,variable.name="strain", value.name="direction")
genes.melt[is.na(genes.melt)] <- 0
genes.melt$strain <- factor(genes.melt$strain,levels=strain.order)

# making a matrix
m <- tidyr::pivot_wider(genes.melt,names_from="gene",values_from="direction")
m <- as.matrix(m[,-1])

# clustering matrix to determine order of genes on x-axis of geom_tile plot
clust <- hclust(dist(t(m)))

# count of all genes
total <- length(unique(genes.melt$gene)); total

# count of shared up/down regulated genes
tmp <- genes.2
tmp[is.na(tmp)] <- 0
tmp$count.up <- apply(tmp, 1, function(x) length(which(x=="2")))
tmp$count.down <- apply(tmp, 1, function(x) length(which(x=="1")))
n.up <- (nrow(tmp[tmp$count.up >= 5,])/total)*100; n.up
n.down <- (nrow(tmp[tmp$count.down >= 5,])/total)*100; n.down

# plotting matrix (Figure 3D)
matrix.plot <- ggplot(genes.melt,aes(gene,strain)) + 
  geom_tile(aes(fill=as.factor(direction)),width=0.5,height=0.5,linewidth=1) + 
  scale_fill_manual(values=c("white",my_colors),labels=c("NS","downregulated","upregulated"),
                    breaks=c("0","1","2"),name="Differential expression in biofilms")+
  scale_x_discrete(limits=colnames(m)[clust$order])+xlab(paste("DE Genes ","(",total," genes total)",sep=""))+
  theme(axis.text.x = element_blank(), axis.ticks= element_blank(), legend.position = "top",
        axis.text.y = element_text(size=12,face="bold"), legend.text = element_text(size=12),
        axis.title.x = element_text(size=12,face="bold"))+
  ylab(NULL) + scale_y_discrete(labels=strain.names, breaks=strain.order)
matrix.plot

ExportPlot(matrix.plot,"../NewFigures/Figure3D",width=12,height=6)

# df of total number of DEGs (split by up/down regulated) per strain
count.df <- data.frame(strain=c("31","31","55","55","345","345","72","72","49","49","540","540"),
                       countType=rep(c("Down","Up"),6),
                       count=c(nrow(data[data$strain == "31" & data$direction == "DOWN",]),
                               nrow(data[data$strain == "31" & data$direction == "UP",]),
                               nrow(data[data$strain == "55" & data$direction == "DOWN",]),
                               nrow(data[data$strain == "55" & data$direction == "UP",]),
                               nrow(data[data$strain == "345" & data$direction == "DOWN",]),
                               nrow(data[data$strain == "345" & data$direction == "UP",]),
                               nrow(data[data$strain == "72" & data$direction == "DOWN",]),
                               nrow(data[data$strain == "72" & data$direction == "UP",]),
                               nrow(data[data$strain == "49" & data$direction == "DOWN",]),
                               nrow(data[data$strain == "49" & data$direction == "UP",]),
                               nrow(data[data$strain == "540" & data$direction == "DOWN",]),
                               nrow(data[data$strain == "540" & data$direction == "UP",])))

# putting strains in correct order
count.df$strain <- factor(count.df$strain,levels=strains)
 
# horizontal bar plot of # of DEGs per strain - to be paired with tree (Figure 2)
count.plot <- ggplot(count.df,(aes(x=strain,y=count))) + geom_col(aes(fill=countType),width=0.5,alpha=0.9)+
  theme_bw() +coord_flip()+theme(axis.text.y=element_blank(),axis.title.y = element_blank(),legend.position = "right",
                                 axis.text.x =element_text(size=14),axis.title.x=element_text(size=14,face="bold"),
                                 legend.text = element_text(size=14),legend.title=element_text(size=14,face="bold"))+
  ylab("Number of DEGs")+ scale_fill_manual(name="\nDifferential expression\nin biofilms",values=my_colors,labels=c("downregulated","upregulated"),breaks=c("Down","Up"))

count.plot

ExportPlot(count.plot,"../NewFigures/Figure2R",width=6,height=4)

# l2fc values by feature type (Figure 8A)
ncrna <- read.delim("Metadata/all_ncRNA.txt",header=F)$V1
sorfs <- read.delim("Metadata/all_sORFs.txt",header=F)$V1
ABvP.result <- data.frame(ABvP.result)
ABvP.result$gene <- rownames(ABvP.result)

ABvP.result <- ABvP.result %>% mutate(type=case_when(gene %in% ncrna ~ "ncRNA",
                                           gene %in% sorfs ~ "sORF",
                                           TRUE ~ "ORF"))
my_colors <- magma(n=4, end = 0.8, begin=0.25)
feature_colors <- c(my_colors[c(2)],"forestgreen",my_colors[c(4)])

plot <- ggplot(ABvP.result,aes(x=type,y=log2FoldChange)) + geom_hline(yintercept=0,linetype="dashed")+
  geom_boxplot(aes(fill=type),alpha=0.8)+
  theme_minimal() + xlab(NULL)+ ylab("Log2 Fold Change")+
  scale_fill_manual(name=NULL,values=feature_colors,breaks=c("ncRNA","ORF","sORF"))+
  theme(axis.text=element_text(size=10,face="bold"),
        axis.text.x = element_blank(),
        axis.title = element_text(size=12),
        strip.text = element_text(size=12,face="bold"))+
  scale_y_continuous(limits=c(-12.5,6),breaks=c(-10,-5,0,5))
plot

stat <- compare_means(log2FoldChange ~ type, data = ABvP.result,p.adjust.method = "BH")
stat <- stat %>% mutate(y.position=c(5,NA,6))
box.stats <- plot + stat_pvalue_manual(stat, label="p.signif",label.size=4)
box.stats

ExportPlot(box.stats,"../NewFigures/Figure8A",width=6,height=6)
