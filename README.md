# Bulk-RNAseq-Fungi
Transcriptomics analysis to important genes and pathways with clinical implications
Ref: https://f1000research.com/articles/5-1408

Author: Ling Wang

```{r eval = FALSE}
library(limma) 
library(edgeR) 
library(tidyr)
library(gplots) 
library(Glimma)
```

# Part 1: Data re-formatting 

```{r eval = FALSE}
RawCount <- read.csv("ECAK_featurecounts.csv")
  head(RawCount)
  dim(RawCount) # 6263 x 18
```

### Reformat the count table: split GeneID and GenBank
```{r eval = FALSE}
RawCount <- RawCount %>%
  separate(Geneid,
           into = c("GeneID", "Genbank"), 
           sep = "(?<=[,])(?=[,Genbank:])" )
  
  head(RawCount)
  dim(RawCount) # 6263 x 18+1
```
  
### some of them don't have "GenBank" but that's okay
```{r eval = FALSE}
RawCount[c(20:23,65:70),1:4]
which(is.na(RawCount)) # See who's got NA
which(is.na(RawCount$Genbank))
```

### Remove "," and ":"
```{r eval = FALSE}
RawCount$GeneID <- gsub(",","",RawCount$GeneID)
RawCount$GeneID <- gsub("GeneID:","",RawCount$GeneID)
RawCount$Genbank <- gsub("Genbank:","",RawCount$Genbank)
  RawCount[c(20:23,65:70),1:4]
  head(RawCount)
  dim(RawCount) # 6263 x 18+1

#write.csv(RawCount,"ECAK_featurecounts_clean.csv", row.names = FALSE)
```

### Read data into DGEList

```{r eval = FALSE}
Count <- RawCount[,c(8:19)] 
  head(Count)
  
Table <- read.csv("EC-AK_SampleTable.csv")
  Table

Gene <- data.frame(Gene=RawCount$GeneID)
  head(Gene)
  dim(Gene) # 6263 x 1
  
x <- DGEList(counts=Count, 
             samples=Table, 
             genes=Gene)
x
```

# Part 2: Organize sample-level information 

```{r eval = FALSE}
SampleName <- c("EC-AK-01","EC-AK-02", "EC-AK-03", 
                "EC-AK-04","EC-AK-05", "EC-AK-06",
                "EC-AK-07","EC-AK-08", "EC-AK-09", 
                "EC-AK-10","EC-AK-11", "EC-AK-12")
colnames(x) <- SampleName

# SampleName
SampleName <- as.factor(c("EC-AK-01","EC-AK-02", "EC-AK-03", 
                          "EC-AK-04","EC-AK-05", "EC-AK-06",
                           "EC-AK-07","EC-AK-08", "EC-AK-09", 
                           "EC-AK-10","EC-AK-11", "EC-AK-12"))
x$samples$SampleName <- SampleName

# SampleType 
SampleType <- as.factor(rep(c("yeast","hyphae"), c(6,6)))
x$samples$SampleType <- SampleType

# Treatment
Treatment <- as.factor(rep(c("ctrl","PYY","ctrl","PYY"), c(3,3,3,3)))
x$samples$Treatment <- Treatment
x$samples

# TotalDescription
TotalDescription <- as.factor(rep(c("yeast_ctrl","yeast_PYY","hyphae_ctrl","hyphae_PYY"), c(3,3,3,3)))
x$samples$TotalDescription <- TotalDescription
x$samples
```


# Prat 3: Data pre-processing 

```{r eval = FALSE}
Before_Preprocessing <- cbind(x$genes,x$counts[,c(1:12)])
#write.csv(Before_Preprocessing,"ECAK-Before_Preprocessing.csv", row.names = FALSE)

#####
##### 1. Remove lowly expressed genes 
#####

# See how many gene counts are ZERO across samples: 33
table(rowSums(x$counts==0)==12)
#FALSE  TRUE 
#6230    33

# Remove them and keep worthwhile genes
keep.exprs <- filterByExpr(x, group=TotalDescription)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
  dim(x) # 6103   12

# What's the default cutoff for filterByExpr? 
# It keeps genes with > 10 read counts 
L <- mean(x$samples$lib.size) * 1e-6
M <- median(x$samples$lib.size) * 1e-6
c(L, M) # 62.60036 62.52972

lcpm.cutoff <- log2(10/M + 2/L)
lcpm.cutoff
# lcpm.cutoff= -2.381771
```
```{r eval = FALSE}
#####
##### 2. Normalize gene expression distributions 
#####

# trimmed mean of M-values (TMM) method in edgeR
# normalisation factors calculated here are used as 
# a scaling factor for the library sizes

x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors
# 0.8094869 0.8250225 0.8131224 0.7713679 0.7479910 0.8001237 1.2844341 1.2702186 1.2646529
# 1.2539405 1.2396073 1.2437430


After_Preprocessing <- cbind(x$genes,x$counts[,c(1:12)])
write.csv(After_Preprocessing,"ECAK-After_Preprocessing.csv", row.names = FALSE)
# reduced #s of genes from 6230+33 to 6103
```

```{r eval = FALSE}
#####
##### 3. Transformations to log scale
#####

cpm <- edgeR::cpm(x)
  head(x$counts)
  head(cpm)
# cpm = count/each library size * 10^6
x$samples$lib.size

#IF getting this message"unable to find an inherited method for function ‘cpm’ for signature ‘"DGEList"’"
#Type the following:
#cpm
#edgeR::cpm 
  
######################################
lcpm <- edgeR::cpm(x, log=TRUE)
  head(lcpm)
  summary(lcpm)
# lcpm = log2(cpm + 2/L) *2 is prior count to avoid log zero
```
  
# Part 4: VennDiagramg 

```{r eval = FALSE}
# See the "common" genes by looking at the filtered count profiles
x$counts[c(1:50),c(1:8)]

par(mfrow=c(2,2))
vennDiagram(x$counts[,c(1:3)], circle.col=c("lightblue1", "lavenderblush3", "darkmagenta"))
title(main="yeast_ctrl",cex.main = 2)

vennDiagram(x$counts[,c(4:6)], circle.col=c("lightblue1", "lavenderblush3", "darkmagenta"))
title(main="yeast_PYY",cex.main = 2)

vennDiagram(x$counts[,c(7:9)], circle.col=c("lightblue1", "lavenderblush3", "darkmagenta"))
title(main="hyphae_ctrl", cex.main = 2)

vennDiagram(x$counts[,c(10:12)], circle.col=c("lightblue1", "lavenderblush3", "darkmagenta"))
title(main="hyphae_PYY",cex.main = 2)


par(mfrow=c(4,3))
```

### See what are those unique genes

```{r eval = FALSE}
Map <- read.csv("ProteinTable21_294796.csv")
Map <- Map[,c(6,7,8,9,11)]
head(Map)
dim(Map)
```


```{r eval = FALSE}
de.common <- which(x$counts[,7]>0 & x$counts[,1]<1)
length(de.common) # 34 unique
unique_hyphae<-data.frame(Gene=head(x$genes$Gene[de.common], n=100))

de.common <- which(x$counts[,10]>0 & x$counts[,4]<1)
length(de.common) # 34 unique
unique_hyphae_pyy<-data.frame(Gene=head(x$genes$Gene[de.common], n=100))

unique_hyphae_ID <- merge(unique_hyphae, Map, by.x = "Gene", by.y = "GeneID")  
  dim(unique_hyphae_ID) #20    10 don't match
  unique_hyphae_ID$Locus 

unique_hyphae_pyy_ID <- merge(unique_hyphae, Map, by.x = "Gene", by.y = "GeneID")  
  dim(unique_hyphae_pyy_ID) # 9 don't match
  unique_hyphae_pyy_ID$Locus 

```


# Part 5: Unsupervised clustering of samples

```{r eval = FALSE}
lcpm <- edgeR::cpm(x, log=TRUE)
head(lcpm)

par(mfrow=c(1,1))
#col.group <- SampleName
#levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
#col.group <- as.character(col.group)

plotMDS(lcpm, labels=SampleName,  
        cex = 1.5, pch=19,
        #col=c(rep("blue4",4), rep("maroon3",4))  )
        col=c("blue4","maroon3","blue4","maroon3",
              "blue4","maroon3","blue4","maroon3")  )
title(main="All samples")


library(Glimma)
glMDSPlot(lcpm, labels=paste(SampleName, SampleName, sep="_"), 
          groups=x$samples[,c(4,5,6)], launch=T)

```


# Part 6: Differential expression analysis 

```{r eval = FALSE}
#####
##### 1. Create a design matrix and contrasts
#####

design <- model.matrix(~0+TotalDescription)
colnames(design) <- gsub("TotalDescription", "", colnames(design))
design

contr.matrix <- makeContrasts(
  yeast_PYYvsyeast_ctrl = yeast_PYY-yeast_ctrl, 
  hyphae_PYYvshyphae_ctrl = hyphae_PYY-hyphae_ctrl,
  hyphae_PYYvsyeast_PYY = hyphae_PYY-yeast_PYY,
  hyphae_ctrlvsyeast_ctrl= hyphae_ctrl-yeast_ctrl,
  levels = colnames(design))
contr.matrix
```

```{r eval = FALSE}
#####
##### 2. Remove heteroscedascity from count data
#####

par(mfrow=c(1,2))
v <- voom(x, design, plot=T)
v

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)

efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")

#summary(decideTests(efit)) # Defalut = p<0.05
```

```{r eval = FALSE}
#####
##### 3. The TREAT method
#####

# Testing significance relative to a fold-change threshold 
# TREAT computes empirical Bayes moderated-t p-values relative to a minimum
# required fold-change threshold.
#  Genes will need to exceed this threshold by some way before being declared 
# statistically significant. 
# It is better to interpret the threshold as
# “the fold-change below which we are definitely not interested in the gene" 
# rather than “the fold-change above which we are interested in the gene"
# https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
# Results: Genes that are significantly DE and have a FC "substaintially" above 
# 1.3 (or lfc=0.3785116) at an FDR (Benjamini Hochberg method) of less than 5% and adj.p.value < 0.05. 


tfit <- treat(vfit, lfc=log2(1.3)) 
dt <- decideTests(tfit) # default: p.value=0.05 + lfc=0 + adjust.method="BH"
summary(dt)  

#hyphae_PYYvshyphae_ctrl
#Down       2
#NotSig   6093
#Up          8
```

```{r eval = FALSE}
#####
##### 4. Plot DE gene distribution
#####
par(mfrow=c(1,2))
plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1],xlim=c(-5,20),ylim=c(-3,4))
plotMD(tfit, column=2, status=dt[,2], main=colnames(tfit)[2],xlim=c(-5,20),ylim=c(-3,4))


# Interactive plots
glMDPlot(tfit, coef=1, status=dt, main=colnames(tfit)[1],
         side.main="ensembl_gene_id", counts=lcpm, groups=SampleName, launch=TRUE)
glMDPlot(tfit, coef=2, status=dt, main=colnames(tfit)[2],
         side.main="ensembl_gene_id", counts=lcpm, groups=SampleName, launch=TRUE)
```

```{r eval = FALSE}
# What are those 32 DE genes in yeast
head(tfit$genes$Gene[which(dt[,1]!=0)], n=1000) 
head(tfit$genes$Gene[which(dt[,1]>0)], n=1000) # Up: 8
head(tfit$genes$Gene[which(dt[,1]<0)], n=1000) # Down: 65

# What are those 60 DE genes in hyphae
head(tfit$genes$Gene[which(dt[,2]!=0)], n=1000) 
head(tfit$genes$Gene[which(dt[,2]>0)], n=1000) # Up: 159
head(tfit$genes$Gene[which(dt[,2]<0)], n=1000) # Down:    118

# What are those 60 DE genes in hyphae-yeast (no PYY)
tfit <- treat(vfit, lfc=log2(20)) 
dt <- decideTests(tfit) # default: p.value=0.05 + lfc=0 + adjust.method="BH"
summary(dt)  

plotMD(tfit, column=4, status=dt[,4], main=colnames(tfit)[4],xlim=c(-5,20),ylim=c(-10,10))

glMDPlot(tfit, coef=4, status=dt, main=colnames(tfit)[4],
         side.main="ensembl_gene_id", counts=lcpm, groups=SampleName, launch=TRUE)


head(tfit$genes$Gene[which(dt[,4]!=0)], n=1000) 
head(tfit$genes$Gene[which(dt[,4]>0)], n=1000) # Up: 13
head(tfit$genes$Gene[which(dt[,4]<0)], n=1000) # Down:    32
```

```{r eval = FALSE}
# hy-yt 
topTreat_HY <- topTreat(tfit,coef=4,n=Inf,p.value=5e-2) 
  dim(topTreat_HY) # 73 X 6
  head(topTreat_HY) 
topTreat_HY_up <- topTreat_HY[topTreat_HY$logFC >= 0, ]
  dim(topTreat_HY_up) # 8
topTreat_HY_down <- topTreat_HY[topTreat_HY$logFC < 0, ]
  dim(topTreat_HY_down) # 65

# Yeast
topTreat_HY_ID <- merge(topTreat_HY, Map, by.x = "Gene", by.y = "GeneID")  
  dim(topTreat_HY_ID) # 72     1 mismatch
  head(topTreat_HY_ID)


write.csv(topTreat_HY_ID,"ECAK_lfc_hyphae_ctrl_vs_yeast_ctrl.csv", row.names = FALSE)
```
