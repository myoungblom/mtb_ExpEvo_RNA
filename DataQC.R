require(DESeq2)
require(dplyr)
require(ggplot2)
require(RColorBrewer)
require(pheatmap)
require(viridis)
require(reshape2)
require(samr)
require(ggpubr)
require(stringr)

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
sampleData <- read.delim("Metadata/mtb_ExpEvo_RNA_metadata.txt", sep="", header = TRUE)

# reformat metadata
sampleTable <- data.frame(SampleName = sampleData$LibraryName,
                          CountsFile = sampleData$CountsFile,
                          Strain = factor(sampleData$Strain),
                          Genotype = factor(sampleData$Genotype),
                          Condition = factor(sampleData$Condition,levels = c("Planktonic","Biofilm")),
                          WetWeight = as.numeric(sampleData$WetWeight),
                          Clade = factor(sampleData$Clade),
                          SampleID = sampleData$SampleID,
                          LibraryPrepBatch = factor(sampleData$LibraryPrepBatch),
                          ReadCount = sampleData$ReadCount,
                          Resequenced = factor(sampleData$Resequenced))

# DE calculations by condition (biofilm vs plantkonic)
Table <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = ".", design = ~ Condition+Genotype+Strain)

# filter out genes with 0 counts in all samples
Table <- Table[rowSums(counts(Table)) > 1,]

# vst transformation
vst <- vst(Table, blind=F)

# PCA plots to look for batch effects (Figure S1)
lib.pca <- DESeq2::plotPCA(vst, intgroup=c("LibraryPrepBatch")) +
  theme_minimal()+ scale_color_manual(name="Library Prep Batch",values=viridis(n=5),
                                      breaks=c("1","2","3","4","5")) +
  theme(axis.text = element_text(size=10),axis.title = element_text(size=12), legend.text = element_text(size=12),
        legend.title = element_text(size=12),legend.position="top")
lib.pca

reseq.pca <- DESeq2::plotPCA(vst, intgroup=c("Resequenced")) +
  theme_minimal()+ scale_color_manual(name="Resequenced?",values=viridis(n=2,end=0.8),
                                      breaks=c("Y","N")) +
  theme(axis.text = element_text(size=10),axis.title = element_text(size=12), legend.text = element_text(size=12),
        legend.title = element_text(size=12),legend.position="top")
reseq.pca

count.pca.data <- DESeq2::plotPCA(vst, intgroup=c("ReadCount","Condition","Genotype"),returnData=T)
count.pca <- ggplot(count.pca.data, aes(x=PC1,y=PC2,color=ReadCount,shape=Condition)) + geom_point(size=4) +
  theme_minimal()+ scale_shape_manual(name="Growth Condition",values=c(16,17),breaks=c("Biofilm","Planktonic"))+
  scale_color_continuous(type = "viridis") +
  theme(axis.text = element_text(size=10),axis.title = element_text(size=12), legend.text = element_text(size=12),
        legend.title = element_text(size=12),legend.position="top")
count.pca

count.pca.2 <- ggplot(count.pca.data, aes(x=PC1,y=PC2,color=ReadCount,shape=Genotype)) + geom_point(size=4) +
  theme_minimal()+ scale_shape_manual(name="Genotype",values=c(16,17),breaks=c("Ancestral","Evolved"))+
  scale_color_continuous(type = "viridis") +
  theme(axis.text = element_text(size=10),axis.title = element_text(size=12), legend.text = element_text(size=12),
        legend.title = element_text(size=12),legend.position="top")
count.pca.2

# calculate distances between samples
dists <- dist(t(assay(vst)))

# matrix of sample distances
distMatrix <- as.matrix(dists)

# convert to data frame, remove 0 distances
dist.df <- as.data.frame(as.table(distMatrix))
colnames(dist.df) <- c("S1","S2","Dist")
dist.df <- dist.df[!(dist.df$Dist == 0),]

# add metadata to data frame
strains1 <- c()
strains2 <- c()
clades1 <- c()
clades2 <- c()
geno1 <- c()
geno2 <- c()
cond1 <- c()
cond2 <- c()

strain.df <- sampleTable %>% mutate(across(everything(), as.character))

for(row in 1:nrow(dist.df)){
  id1 <- dist.df$S1[row]
  id2 <- dist.df$S2[row]
  strains1 <- c(strains1,strain.df$Strain[strain.df$SampleName == id1])
  strains2 <- c(strains2,strain.df$Strain[strain.df$SampleName == id2])
  clades1 <- c(clades1,strain.df$Clade[strain.df$SampleName == id1])
  clades2 <- c(clades2,strain.df$Clade[strain.df$SampleName == id2])
  geno1 <- c(geno1,strain.df$Genotype[strain.df$SampleName == id1])
  geno2 <- c(geno2,strain.df$Genotype[strain.df$SampleName == id2])
  cond1 <- c(cond1,strain.df$Condition[strain.df$SampleName == id1])
  cond2 <- c(cond2,strain.df$Condition[strain.df$SampleName == id2])
}

dist.df$S1_strain <- strains1
dist.df$S2_strain <- strains2
dist.df$S1_clade <- clades1
dist.df$S2_clade <- clades2
dist.df$S1_geno <- geno1
dist.df$S2_geno <- geno2
dist.df$S1_cond <- cond1
dist.df$S2_cond <- cond2

dist.df <- dist.df %>% 
  mutate(Strain = case_when(S1_strain == S2_strain ~ S1_strain,
                          S1_strain != S2_strain ~ "diff")) %>%
  mutate(Clade = case_when(S1_clade == S2_clade ~ S1_clade,
                          S1_clade != S2_clade ~ "diff")) %>%
  mutate(Cond = case_when(S1_cond == S2_cond ~ S1_cond,
                          S1_cond != S2_cond ~ "diff")) %>%
  mutate(Geno = case_when(S1_geno == S2_geno ~ S1_geno,
                          S1_geno != S2_geno ~ "diff"))

dist.df <- dist.df[,c("S1","S2","Dist","Cond","Geno","Strain","Clade")]

dist.df <- dist.df %>% mutate(Rep = case_when(!(Cond == "diff" | Geno == "diff" | Strain == "diff") ~ "YES",
                                              (Cond == "diff" | Geno == "diff" | Strain == "diff") ~ "NO"))

dist.df$Dist <- as.numeric(dist.df$Dist)

dist.df<- dist.df[!duplicated(apply(dist.df,1,function(x) paste(sort(x),collapse=''))),]

# only within condition & genotype comparisons
dist.df.wi <- dist.df[dist.df$Cond != "diff" & dist.df$Geno != "diff",]

# plot distances of biological replicates (Figure S2)
rep.dist <- ggplot(dist.df.wi) + geom_boxplot(aes(x=Rep,y=Dist,fill=Rep)) + theme_bw()+
  xlab("Replicates?") + ylab("Pairwise Distance") + theme(legend.position = "none")
rep.dist

rep.stats <- compare_means(data=dist.df.wi, formula = Dist ~ Rep, p.adjust.method = "bonferroni")
rep.stats

# remove distances between replicates
dist.df.wi.norep <- dist.df.wi[dist.df.wi$Rep == "NO",]

# plot sample distances by genotype (Figure S5B)
geno.dist <- ggplot(dist.df.wi.norep) + geom_boxplot(aes(x=Cond,y=Dist,fill=Geno)) + theme_bw() +
  xlab("Condition") + ylab("Pairwise Distance") + 
  scale_fill_manual(name="Genotype",values=c("#7CAE00","#C77CFF"))
geno.dist

stats <- compare_means(data=dist.df.wi.norep, formula = Dist ~ Geno, group.by = "Cond", p.adjust.method = "bonferroni")
stats

                              