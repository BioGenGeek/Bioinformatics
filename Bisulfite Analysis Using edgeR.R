### Bisulfite sequencing data analysis 
### Identifying differential methylated CpG Loci
library(edgeR)
library(GEOquery)
library(tidyverse)
library(data.table)

### Bisulfite sequencing data: GSE86297
gse<- getGEO('GSE86297', GSEMatrix = TRUE, destdir = 'E:/Codes')

targets <- read.delim("targets.txt", stringsAsFactors=FALSE)
Sample1 <- read.delim(
  file="E:/GSE86297_RAW/GSM2299710_RRBS_40-45oocyte_LibA.cov.txt.gz",
  header=FALSE)
setnames(Sample1, c("Chr", "start", "end", "methylation%", 
                    "methylated", "un_methylated"))
**load all six samples as a DGEList**
#yall <- readBismark2DGE(targets$File, sample.names=targets$Sample)
smpl<-c('GSM2299710_RRBS_40-45oocyte_LibA.cov.txt.gz',
        'GSM2299711_RRBS_40-45oocyte_LibB.cov.txt.gz',
        'GSM2299712_RRBS_50-55oocyte_LibA.cov.txt.gz',
        'GSM2299713_RRBS_50-55oocyte_LibB.cov.txt.gz',
        'GSM2299714_RRBS_60-65oocyte_LibA.cov.txt.gz',
        'GSM2299715_RRBS_60-65oocyte_LibB.cov.txt.gz')
sample_names<-c( '40-45-A', '40-45-B','50-55-A','50-55-B','60-65-A','60-65-B')
setwd('E:/GSE86297_RAW')
alldata<- readBismark2DGE(smpl, sample.names=sample_names)
dim(alldata)
counts <- alldata$counts
samples <- alldata$samples
genes <- alldata$genes
table(alldata$genes$Chr)

### remove the mitochondrial genes
alldata <- alldata[alldata$genes$Chr!="MT", ]
table(alldata$genes$Chr)

### sort the DGEList for all loci in genomic order, from Chr 1 to Chr Y
ChrNames <- c(1:19, "X", "Y")
alldata$genes$Chr <- factor(alldata$genes$Chr, levels = ChrNames)
chr_order <- order(alldata$genes$Chr, alldata$genes$Locus)
alldata <- alldata[chr_order,]
table(alldata$genes$Chr)

### annotate the CpG loci with the identity of the nearest gene
# Search for the gene transcriptional start site (TSS) closest to each CpGs
TSS <- nearestTSS(alldata$genes$Chr, alldata$genes$Locus, species="Mm")

# add TSS information into DGEList$genes
head(alldata$genes)
alldata$genes$EntrezID <- TSS$gene_id
alldata$genes$Symbol <- TSS$symbol
alldata$genes$Width <- TSS$width
alldata$genes$Strand <- TSS$strand
alldata$genes$Distance <- TSS$distance # genomic distance from the CpG to the TSS
head(alldata$genes)

### Filtering and normalization
# get the total read coverage at each CpG site for each sample
Cytosine <- gl(2, 1, ncol(alldata), labels=c("Me", "Un"))
Me <- alldata$counts[, Cytosine=="Me"]
Un <- alldata$counts[, Cytosine=="Un"]
Coverage <- Me + Un # sum up the counts of methylated and unmethylated reads
head(Coverage)

# Filtering: remove CpG loci that have low coverage
# A CpG site to have a total count (me & un) of at least 8 in every sample
HasCoverage <- rowSums(Coverage >= 8) == 6
table(HasCoverage)
HasBoth <- rowSums(Me) > 0 & rowSums(Un) > 0
table(HasBoth)
table(HasCoverage, HasBoth)

# filter the DGEList
y <- alldata[HasCoverage & HasBoth, , keep.lib.sizes=FALSE]

### No normalization: set the library sizes equal for each pair of libraries
y$samples
TotalLibSize <- y$samples$lib.size[Cytosine=="Me"] +
  y$samples$lib.size[Cytosine=="Un"]
TotalLibSize <- TotalLibSize/2
y$samples$lib.size <- rep(TotalLibSize, each=2)
y$samples

### Data exploration: plot methylation level (M-value) of the CpG sites
Me <- y$counts[, Cytosine=="Me"]
Un <- y$counts[, Cytosine=="Un"]
M <- log2(Me + 2) - log2(Un + 2) # +2 to avoid logarithms of zero
#colnames(M) <- targets$Sample
colnames(M) <- sample_names
plotMDS(M, col=rep(1:3, each=2), main="M-values")

### Design matrix
designSL <- model.matrix(~0+Group, data=target)
designSL

 # expand to the full design matrix modeling the sample and methylation effects
design <- modelMatrixMeth(designSL) 

design # six columns represent the sample coverage effects, 
# last three columns represent the methylation levels in the three groups
### Dispersion estimation
# CpG methylation in chromosome 1, subset chromosome 1
y1 <- y[y$genes$Chr==1, ]

# estimate the NB dispersion for each CpG site
y1 <- estimateDisp(y1, design=design, trend="none") # not consider trend 
summary(y1$prior.df) 

# all the CpG-wise dispersions are exactly the same as the common dispersion
### Differential methylation analysis at CpG loci between different groups
fit <- glmFit(y1, design)
contr <- makeContrasts(Group60vs40 = Group60um - Group40um, levels=design)
lrt <- glmLRT(fit, contrast=contr)
topTags(lrt)
summary(decideTests(lrt))
plotMD(lrt)
