require(DESeq2)
require(dplyr)
require(ggplot2)
require(RColorBrewer)
require(pheatmap)
require(viridis)
require(reshape2)

#####
# Analysis of M. tb RNA-seq data from evolved populations, using DESeq2 and counts from HTSeq.
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

# subset to only evolved samples
sampleTable <- sampleTable[sampleTable$Genotype == "Evolved",]

# DE calculations by condition (biofilm vs plantkonic)
DESeq2Table <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = ".", design = ~ Condition)

# filter out genes with 0 counts in all samples
DESeq2Table <- DESeq2Table[rowSums(counts(DESeq2Table)) > 1,]

# run DESeq analysis
DESeq2Table <- DESeq(DESeq2Table)
resultsNames(DESeq2Table)

# r-log transformation
rld <- rlogTransformation(DESeq2Table, blind=FALSE)

# pca plot colored by growth condition (Figure 4B)
evo.pca.plot <- DESeq2::plotPCA(rld, intgroup=c("Condition")) +
  theme_minimal()+ scale_color_manual(name="Growth Condition",values=c("#150E37FF","#D1426FFF"),
                                      breaks=c("Biofilm","Planktonic"))
evo.pca.plot

# results of biofilm vs planktonic comparison
result <- results(DESeq2Table, alpha = 0.05,lfcThreshold = 0, contrast=c("Condition","Biofilm","Planktonic"))
# number of significantly DE genes
sum(result$padj <0.05, na.rm=TRUE)
# subset all genes for only significant results (Supplementary Data 2 - AllPopulations_EBvP)
subset <- subset(result,padj<0.05)
write.csv("../2023.02.28_newFigs/DEFiles/all_EBvP_DE.csv",x=subset)

# heatmap of all genes (Figure S2B)
my_colors <- magma(n=4, end = 0.8, begin=0.1)
condition <- my_colors[c(1,3)]
genotype <- my_colors[c(2,4)]
names(condition) <- c("Biofilm","Planktonic")
anno_colors <- list(Condition=condition)
BvPCounts <- assay(rld)[rownames(result),]
BvPCounts <- BvPCounts - rowMeans(BvPCounts)
anno <- as.data.frame(colData(rld)[c("Condition")])
custom.palette <- c("red","white","blue")
colors <- colorRampPalette(rev(custom.palette))(20)

de.heatmap <- pheatmap(BvPCounts, color = colors, show_rownames = F, annotation_col = anno,show_colnames = F, 
                       cutree_cols = 2,treeheight_row = 0, annotation_colors = anno_colors,
                       kmeans_k = NA, border_color = NA, scale= "column", fontsize = 12)

de.heatmap

#volcano plot (Figure 4C)
voldata <- data.frame(result)
voldata <- voldata %>% mutate(sig=case_when(padj < 0.05 & log2FoldChange > 0 ~ "Significant-Up",
                                            padj < 0.05 & log2FoldChange < 0 ~ "Significant-Down",
                                            padj > 0.05 ~ "Not significant"))

custom.palette <- colorRampPalette(rev(c("red","white","blue")))(20)
my_colors <- c("grey",custom.palette[15],custom.palette[5])

vol <- ggplot(voldata,aes(x=log2FoldChange,y=-log10(padj))) + geom_point(aes(color=sig),alpha=0.8) + 
  theme_minimal()+xlim(-11,11)+ylim(0,80)+
  scale_color_manual(name=NULL,values=my_colors,breaks=c("Not significant","Significant-Up","Significant-Down"),
                     labels=c("Not significant","Upregulated","Downregulated"))+
  theme(legend.position = "top")
vol
ExportPlot(vol,"../2023.02.28_newFigs/IndvFigures/all_EBvP_volanco",width=6,height=6)

# separate by individual strains (Supplementary Data 2 - EBvP for each strain)
strains <- c("31","49","55","72","345","540")

for (strain in strains){
  st <- sampleTable[sampleTable$Strain == strain,]
  t <- DESeqDataSetFromHTSeqCount(sampleTable = st, directory = ".", design = ~ Condition)
  t <- t[rowSums(counts(t)) > 1,]
  t <- DESeq(t)
  r <- results(t, alpha = 0.05, lfcThreshold = 0, contrast=c("Condition","Biofilm","Planktonic"))
  sub <- subset(r, padj < 0.05)
  write.csv(sub,paste("../2023.02.28_newFigs/DEFiles/",paste(strain,"evoBvP",sep="_"),".csv",sep=""))
}

# reading individual data files
x31 <- data.frame(read.csv("DEFiles/EBvP/31_evoBvP.csv",header=T))
x49 <- data.frame(read.csv("DEFiles/EBvP/49_evoBvP.csv",header=T))
x55 <- data.frame(read.csv("DEFiles/EBvP/55_evoBvP.csv",header=T))
x72 <- data.frame(read.csv("DEFiles/EBvP/72_evoBvP.csv",header=T))
x345 <- data.frame(read.csv("DEFiles/EBvP/345_evoBvP.csv",header=T))
x540 <- data.frame(read.csv("DEFiles/EBvP/540_evoBvP.csv",header=T))

# reformatting data into one df
data <- rbind(#data.frame(strain="31", gene = x31[x31$log2FoldChange >0,]$X, direction="UP"),
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

# overlap between strains
genes.2[is.na(genes.2)] <- 0
total <- length(unique(genes.melt$gene))
up.counts <- c()
down.counts <- c()
for(x in 1:nrow(genes.2)){
  up.count <- sum(genes.2[x,][-1] == 2)
  up.counts <- c(up.counts,up.count)
  down.count <- sum(genes.2[x,][-1] == 1)
  down.counts <- c(down.counts,down.count)
}

up.shared <- sum(up.counts >= 5)
down.shared <- sum(down.counts >= 5)

print(paste("Number of upregulated genes shared by 5+ strains: ",(up.shared/total)*100,"%"))
print(paste("Number of downregulated genes shared by 5+ strains: ",(down.shared/total)*100,"%"))

# plotting matrix (Figure 4D)
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

ExportPlot(matrix.plot,"../2023.02.28_newFigs/IndvFigures/EBvP_matrix",width=12,height=6)

################################################################
#### FIGURE S3 -- treating MT31 outlier as a biofilm sample ####
################################################################

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

# edit outlier sample condition
sampleTable$Condition[sampleTable$SampleName == "3-E-P"] <- "Biofilm"

# DE calculations by condition (biofilm vs plantkonic)
DESeq2Table <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = ".", design = ~ Condition)

# filter out genes with 0 counts in all samples
DESeq2Table <- DESeq2Table[rowSums(counts(DESeq2Table)) > 1,]

# run DESeq analysis
DESeq2Table <- DESeq(DESeq2Table)
resultsNames(DESeq2Table)

# r-log transformation
rld <- rlogTransformation(DESeq2Table, blind=FALSE)

# pca plot colored by growth condition (Figure S3A)
evo.pca.plot <- DESeq2::plotPCA(rld, intgroup=c("Condition")) +
  theme_minimal()+ scale_color_manual(name="Growth Condition",values=c("#150E37FF","#D1426FFF"),
                                      breaks=c("Biofilm","Planktonic"))
evo.pca.plot

# results of biofilm vs planktonic comparison
result <- results(DESeq2Table, alpha = 0.05,lfcThreshold = 0, contrast=c("Condition","Biofilm","Planktonic"))
# number of significantly DE genes
sum(result$padj <0.05, na.rm=TRUE)
# subset all genes for only significant results
subset <- subset(result,padj<0.05)

# heatmap of all genes (Figure S3C)
my_colors <- magma(n=4, end = 0.8, begin=0.1)
condition <- my_colors[c(1,3)]
genotype <- my_colors[c(2,4)]
names(condition) <- c("Biofilm","Planktonic")
anno_colors <- list(Condition=condition)
BvPCounts <- assay(rld)[rownames(result),]
BvPCounts <- BvPCounts - rowMeans(BvPCounts)
anno <- as.data.frame(colData(rld)[c("Condition")])
custom.palette <- c("red","white","blue")
colors <- colorRampPalette(rev(custom.palette))(20)

de.heatmap <- pheatmap(BvPCounts, color = colors, show_rownames = F, annotation_col = anno,show_colnames = F, 
                       cutree_cols = 2,treeheight_row = 0, annotation_colors = anno_colors,
                       kmeans_k = NA, border_color = NA, scale= "column", fontsize = 12)

de.heatmap

#volcano plot (Figure S3B)
voldata <- data.frame(result)
voldata <- voldata %>% mutate(sig=case_when(padj < 0.05 & log2FoldChange > 0 ~ "Significant-Up",
                                            padj < 0.05 & log2FoldChange < 0 ~ "Significant-Down",
                                            padj > 0.05 ~ "Not significant"))

custom.palette <- colorRampPalette(rev(c("red","white","blue")))(20)
my_colors <- c("grey",custom.palette[15],custom.palette[5])

vol <- ggplot(voldata,aes(x=log2FoldChange,y=-log10(padj))) + geom_point(aes(color=sig),alpha=0.8) + 
  theme_minimal()+xlim(-11,11)+
  scale_color_manual(name=NULL,values=my_colors,breaks=c("Not significant","Significant-Up","Significant-Down"),
                     labels=c("Not significant","Upregulated","Downregulated"))+
  theme(legend.position = "top")
vol

# separate by individual strains
strains <- c("31","49","55","72","345","540")

for (strain in strains){
  st <- sampleTable[sampleTable$Strain == strain,]
  t <- DESeqDataSetFromHTSeqCount(sampleTable = st, directory = ".", design = ~ Condition)
  t <- t[rowSums(counts(t)) > 1,]
  t <- DESeq(t)
  r <- results(t, alpha = 0.05, lfcThreshold = 0, contrast=c("Condition","Biofilm","Planktonic"))
  sub <- subset(r, padj < 0.05)
  write.csv(sub,paste("DEFiles/EBvP_MT31outlier/",paste(strain,"evoBvP",sep="_"),".csv",sep=""))
}

# reading individual data files
x31 <- data.frame(read.csv("DEFiles/EBvP_MT31outlier/31_evoBvP.csv",header=T))
x49 <- data.frame(read.csv("DEFiles/EBvP_MT31outlier/49_evoBvP.csv",header=T))
x55 <- data.frame(read.csv("DEFiles/EBvP_MT31outlier/55_evoBvP.csv",header=T))
x72 <- data.frame(read.csv("DEFiles/EBvP_MT31outlier/72_evoBvP.csv",header=T))
x345 <- data.frame(read.csv("DEFiles/EBvP_MT31outlier/345_evoBvP.csv",header=T))
x540 <- data.frame(read.csv("DEFiles/EBvP_MT31outlier/540_evoBvP.csv",header=T))

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

# overlap between strains
genes.2[is.na(genes.2)] <- 0
total <- length(unique(genes.melt$gene))
up.counts <- c()
down.counts <- c()
for(x in 1:nrow(genes.2)){
  up.count <- sum(genes.2[x,][-1] == 2)
  up.counts <- c(up.counts,up.count)
  down.count <- sum(genes.2[x,][-1] == 1)
  down.counts <- c(down.counts,down.count)
}

up.shared <- sum(up.counts >= 5)
down.shared <- sum(down.counts >= 5)

print(paste("Number of upregulated genes shared by 5+ strains: ",(up.shared/total)*100,"%"))
print(paste("Number of downregulated genes shared by 5+ strains: ",(down.shared/total)*100,"%"))

# plotting matrix (Figure S3D)
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
