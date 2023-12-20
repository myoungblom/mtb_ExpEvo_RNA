require(DESeq2)
require(dplyr)
require(ggplot2)
require(vsn)
require(RColorBrewer)
require(pheatmap)
require(apeglm)
require(genefilter)
require(hexbin)
require(viridis)

#####
# Analysis of M. tb RNA-seq data from planktonic populations, using DESeq2 and counts from HTSeq.
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

# subset to only biofilm samples
P.table <- sampleTable[sampleTable$Condition == "Planktonic",]

# create sample table for DESeq2 input
P.DESeq <- DESeqDataSetFromHTSeqCount(sampleTable = P.table, directory = ".",design = ~ Genotype)

# Planktonic PCA with arrow
P.vst <- vst(P.DESeq, blind=FALSE)
PCAdata <- DESeq2::plotPCA(P.vst, intgroup=c("Clade","Genotype"),returnData=T)
#write.table(PCAdata,file="DEFiles/PEvA/P_EvA_PCAdata.tsv",sep="\t")

# PCA data file edited to connect ancestral and evolved samples
PCAdata2 <- read.table("DEFiles/PEvA/P_EvA_PCAdata_v2.txt",sep="\t",header=T)
rv <- rowVars(assay(P.vst))
select <- order(rv, decreasing = TRUE)[seq_len(min(500,length(rv)))]
pca <- prcomp(t(assay(P.vst)[select, ]))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
PC1var <- round(percentVar[1]*100)
PC2var <- round(percentVar[2]*100)

# make plot with arrows (Figure 5B)
colors <- c(plasma(2,begin=0.4,end = 0.8))
arrows <- arrow(ends="last", length=unit(0.1,"inches"),angle=20,type="closed")

PCAplot <- ggplot()+geom_point(data=PCAdata2,aes(x=PC1,y=PC2,color=Clade,shape=Genotype),size=4,alpha=0.9)+
  geom_curve(aes(x=PC1,y=PC2,xend=PC12,yend=PC22,color=Clade),alpha=0.9,data=PCAdata2,
             arrow=arrows,curvature = 0.1,na.rm = T,show.legend = F,linetype="dashed")+
  theme_minimal()+xlab(paste0("PC1: ",PC1var,"%"))+ylab(paste0("PC2: ",PC2var,"%"))+
  scale_color_manual(name="Sub-lineage",labels=c("L4.9","L4.4"),breaks=c("L4.9","L4.4"),values=colors)+
  scale_shape_manual(name="Genotype",breaks=c("Ancestral","Evolved"), values=c(10,16))+
  theme(legend.position = "top",axis.text=element_text(size=10),axis.title=element_text(size=12),
        legend.text=element_text(size=12))

PCAplot

# separate by individual strains
strains <- c("31","49","55","72","345","540")

for (strain in strains){
  st <- P.table[P.table$Strain == strain,]
  t <- DESeqDataSetFromHTSeqCount(sampleTable = st, directory = ".", design = ~ Genotype)
  t <- t[rowSums(counts(t)) > 1,]
  t <- DESeq(t)
  r <- results(t, alpha = 0.05, lfcThreshold = 0, contrast=c("Genotype","Evolved","Ancestral"))
  #write.csv(r,paste("DEFiles/PEvA/",paste(strain,"PEvA","allGenes",sep="_"),".csv",sep=""),quote=F)
  sub <- subset(r, padj < 0.05)
  #write.csv(sub,paste("DEFiles/PEvA/",paste(strain,"PEvA",sep="_"),".csv",sep=""),quote=F)
}
