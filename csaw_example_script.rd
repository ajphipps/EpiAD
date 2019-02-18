---
title: "H3K27ac_AD_12m_pairwise_analysis_csaw"
author: "Andrew Phipps"
date: "11 July 2018"
output:
  html_document: default
  pdf_document: default
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
Download and install the required packages
```{r, include=FALSE}
#download and install required packages - note you need R 3.3.3 or later
#source("https://bioconductor.org/biocLite.R")
#biocLite("edgeR")
#biocLite("csaw")
#biocLite("DiffBind")
#install.packages("tidyverse")
#install.packages("dplyr")
#biocLite("Gviz")
```
We will now load the required packages into R
```{r}
#load required packages into R
#tutorial link https://www.bioconductor.org/help/course-materials/2015/BioC2015/csaw_lab.html
library(csaw)
library(edgeR)
library(rtracklayer)
library(ggplot2)
library(dplyr)
library(GenomicRanges)

```
Importing the sample files for analysis - note I have included the input control, but not needed for this analysis
For this analysis we will only use a single 

```{r}
#import required bam files - for now stick to single timepoint WT vs TG for K27ac - but not sure if it would be better to load all data
#for multiple marks at a single time point or import all data from the experiment 
#load samples in from marked duplicated files

K27ac_bam.files <- file.path("/NGS_Data/Andrew/Pseudoreplicate_BAM/H3K27ac_pseudoreplicate_bam/", c("H3K27ac_WT12_pseudo1.sorted.bam", "H3K27ac_WT12_pseudo2.sorted.bam", "H3K27ac_TG12_pseudo1.sorted.bam", "H3K27ac_TG12_pseudo2.sorted.bam"))
```
We will now make a design matrix for these samples (keep the same order as samples)
```{r}
#design matrix
genotype <- factor(c('WT12','WT12','TG12','TG12'), levels = c('WT12','TG12'))
design <- model.matrix(~genotype)
colnames(design) <- c("intercept", "genotype")
#making blacklist granges with Rtracklayer
library(rtracklayer)
gr_blacklist = import("/NGS_Data/Andrew/bigwig/mm10.blacklist.bed")
```
Also import mouse mm10 genomic features
```{r}
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
broads <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)

```

This is now where we set up the parameters for importing the data and bring it inot a ranged summarized experiment
```{r}
###############################################################################
#setting up parameters
frag.length <- 250 #this is the fragment length from sonication
window.width <- 150 #will partition into windows of 150bp - for TF set small window, but for histone marks set window relevent \
##to mark - around 150bp windows - min size of a nucleosome? see 2.5 - Choosing appropriate window size
spacing <- 50
#discard mm10 blacklist regions - remember to include encode blacklist regions #probability mapped incorrectly: 0.01 (mapq 20)
parameters <- readParam(minq=30, discard = gr_blacklist)
data <- windowCounts(K27ac_bam.files, ext=frag.length, width =window.width, param = parameters)
```
Visualisation of data to ensure it has loaded in correctly:
```{r}
#can visualise the data with 
rowRanges(data)
```
and visualisation of counts:
```{r}
#visualise counts with 
head(assay(data), n=100)
```
Confirmation of fragment size:
```{r}
#can visualise to pick fragment size with the following plot - however you need to have marked (not removed) duplicates with Picard
#Sharp spike = fragment length. Low curve = potentially poor IP efficiency (low input from my samples)
max.delay <- 500
dedup.on <- reform(parameters, dedup=TRUE)
plot1 <- correlateReads(K27ac_bam.files, max.delay, param = dedup.on)
plot(0:max.delay, plot1, type='l', ylab="CCF", xlab="Delay (bp)")
#can quantitate the fragment size with the following: 
maximizeCcf(plot1)
```
Performing filtering steps - filtering by global enrichment to remove background (see chapter 3)
Note that we are keeping any windows with greater than 3 fold change against background binning
```{r}
#filtering from the 'negative binomial' - log transform of NB is referred to as abundance
##to remove uninteresting regions and lower computational complexity
#a better option for our analysis is to filter via global enrichement
bin.size <- 2000L 
binned1 <- windowCounts(K27ac_bam.files, bin=TRUE, width=bin.size, param=parameters)
filtered.stat <- filterWindows(data, background=binned1, type="global")
#keep samples with a fold change of greater than 4 (3 by default in csaw guide)
filtered.keep <- filtered.stat$filter >log2(4)
```
Total number of filtered samples
```{r}
#sum of filtered samples
sum(filtered.keep)
#make filtered.data array for downstream analysis
filtered.data <- data[filtered.keep,]
```
Generate histogram to show the majority of sites are filtered as background
```{r}
#visualise fold change to confirm that the bulk of background sites are removed by filtering
par(mfrow=c(1,1)) 
hist(filtered.stat$back.abundances, xlab="adjusted bin log-CPM", breaks=100,
     main="", col="grey80", xlim=c(min(filtered.stat$back.abundances), 0))
global.bg <- filtered.stat$abundances -filtered.stat$filter
abline(v=global.bg[1], col="red", lwd=2)
abline(v=global.bg[1]+log2(4), col="blue", lwd=2)
legend("topright", lwd=2, col=c('red', 'blue'), legend=c("Background", "Threshold"))
```
Now use TMM binned normalisation for elimination of composition bias
```{r}
#elimination of composition bias with TMM on binned reads of 10000 bp

binned <- windowCounts(K27ac_bam.files, bin=TRUE, width=10000, param=parameters)
filtered.data <- normOffsets(binned, se.out=filtered.data)
filtered.data$norm.factors
#You can test multiple normalisation windows here: too small = low counts and loss of DB, too large and DB will
##be in the same window as background
#demo <- windowCounts(K27ac_bam.files, bin=TRUE, width=5000, param=parameters)
#normOffsets(demo, se.out=FALSE)
# [1] 0.9748893 1.0295585 0.8987019 1.0386579 1.0815877 0.8709669 0.9466737 1.0718893 1.0981895 1.0167509
#demo <- windowCounts(K27ac_bam.files, bin=TRUE, width=15000, param=parameters)
#normOffsets(demo, se.out=FALSE)
#[1] 0.9847623 1.0302603 0.9183524 1.0549877 1.0909148 0.8883423 0.9719159 1.0686444 1.1166129 0.9051679
```
Generate MA plots to visualise composition bias and normalisation. A vertical shift in the bar might indicate composition bias. We ideally want the composition factors (line) to pass through the centre of the dense cloud
```{r}
#visualisation of normalisation with MA plots - generate for each sample 
#vertical shift in the bars might indicate composition bias - ideally want comp factors (line) to pass through
##centre of the cloud
par(mfrow=c(3,3), mar=c(5,4,2,1.5))
adj.counts <- cpm(asDGEList(binned), log=TRUE)
normfacs <- filtered.data$norm.factors
for (i in seq_len(length(K27ac_bam.files)-1)) {
  cur.x <- adj.counts[,1]
  cur.y <- adj.counts[,1+i]
  smoothScatter(x=(cur.x+cur.y)/2+6*log2(10), y=cur.x-cur.y,
  xlab="A", ylab="M", main=paste("1 vs",i+1))
all.dist <-diff(log2(normfacs[c(i+1, 1)]))
abline(h=all.dist, col="red")
}
```
Now we have set up the samples and everything is ready for testing of differential binding
```{r}
#testing for Diffbinding
#need: filtered.data and filtered.data, original data: K27ac_bam.files, and design matrix
#setting up the data
y<- asDGEList(filtered.data)
#experimental design and setting up data
design
#stabilising estimates with empyrical bayes
y<- estimateDisp(y, design)
summary(y$trended.dispersion)

fit <-glmQLFit(y, design, robust=TRUE)
summary(fit$var.post)
```
We can visualise the empyrical bayes stabilisation for NB dispersion
```{r}
#visualisation of EB stabilisation biological coefficient of variation for NB dispersion - see pg 42

par(mfrow=c(1,2))
o<-order(y$AveLogCPM)
plot(y$AveLogCPM[o],sqrt(y$trended.dispersion[o]), type="l", lwd=2,
                         ylim=c(0,1), xlab=expression("Ave."~Log[2]~"CPM"),
                         ylab=("biological coefficient of variation"))
     plotQLDisp(fit)
#my filtering may not have been strong enough - might need greater fold change to reduce the crap I am seeing here
     #alternatively it might be due to increased variation I see between my samples here? 

summary(fit$df.prior)
```
Now we test for differential binding
#visualise with MDS plots:
#plotMDS()
```{r}
#Testing for differential binding windows
results <- glmQLFTest(fit, contrast=c(0,1))
head(results$table)

#assign p-values to co-ordinates
rowData(filtered.data) <-cbind(rowData(filtered.data),results$table)
```
By utilising MDS plots we can examine similarity between samples for the top 100,5000,1000,5000 called sites
```{r}
#examine replicate similarity with MDS plots
par(mfrow=c(2,2), mar=c(5,4,2,2))
adj.counts<-cpm(y,log=TRUE)
for(top in c(100,500,1000,5000)) {
  out <- plotMDS(adj.counts, main=top, col=c("blue","blue", "red", "red"),
                 labels=c("WT12", "WT12", "TG12", "TG12"), top=top)
}
```
Now we can correct these sites for false discovery rate using Benjamini-Hochbeg method per p-value. This will correct for regions, but is less conservative than Bonferroni.

Clustering with external cues - like mm10 genomic features and background

```{r}
#cluster windows into regions without external information/cues - note max.width to limit size (6.2.2)
#tolerance is the min distance for two binding sites to be treated as separate events
mergedwindowsK27ac <- mergeWindows(rowRanges(filtered.data), tol = 1000, max.width = 8000L)
mergedwindowsK27ac$region
#assigning combined P value for merged windows
p.mergedwindowsK27ac <- combineTests(mergedwindowsK27ac$id, results$table)
#check to see if most clusters are an acceptable size, if there are huge clusters we need to improve our filtering or limit
summary(width(mergedwindowsK27ac$region))
```

Cluster with external cues/prior information - the mm10 DB
```{r}
olap <- findOverlaps(broads, rowRanges(filtered.data))
tabbroad <- combineOverlaps(olap, results$table)
head(tabbroad[!is.na(tabbroad$PValue),],n=100)
```
merged <- mergeWindows(rowRanges(filtered.data), )
```{r}
is.sig.gene <- tabbroad$FDR <= 0.05
table(tabbroad$direction[is.sig.gene])

tabbroad %>% arrange(PValue) %>% head(n = 10)
tabbroad %>% arrange(FDR) %>% head(n = 10)
```

Try quick and dirty clustering to see if we can then annotate genomic features
```{r}
merged <- mergeWindows(rowRanges(filtered.data), tol = 1000L, max.width = 10000L)
merged$region
tabcom <- combineTests(merged$id, results$table)
head(tabcom)
summary(width(merged$region))
tabcom %>% arrange(PValue) %>% head(n = 10)
tabcom %>% arrange(FDR) %>% head(n = 10)
table(tabcom$direction)
#summarising the direction of DB per cluster based on combined p value
is.sig.region <- tabcom$FDR <=0.05
table(tabcom$direction[is.sig.region])
mcols(merged$region) <- tabcom

#tabbest for plotting
tab.best <- getBestTest(merged$id, results$table)
head(tab.best)


```


Gene based annotation for broads dataset
```{R}
#This will ID overlaps between regions and annotated genomic features:
#promoter region will be 3KB upstream and 1Kbp downstream of TSS for that gene. 
#exonic features within dist will also be reported - dist = 5000 bp 
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(GenomicRanges)
library(csaw)
mcols(broads) <- tabbroad
anno.broad <- detailRanges(broads, txdb=TxDb.Mmusculus.UCSC.mm10.knownGene, orgdb=org.Mm.eg.db, promoter=c(3000,1000), dist=5000)
anno.ranges.broad <- detailRanges(broads, txdb=TxDb.Mmusculus.UCSC.mm10.knownGene, orgdb=org.Mm.eg.db)
head(anno.ranges.broad$overlap)

```
gene based annotation for basic clustering
```{r}
#This will ID overlaps between regions and annotated genomic features:
#promoter region will be 3KB upstream and 1Kbp downstream of TSS for that gene. 
#exonic features within dist will also be reported - dist = 5000 bp 
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(GenomicRanges)
library(csaw)

anno.merged <- detailRanges(merged$region, txdb=TxDb.Mmusculus.UCSC.mm10.knownGene, orgdb=org.Mm.eg.db, promoter=c(3000,1000), dist=5000)
head(anno.merged$overlap)
merged$region$overlap <-anno.merged$overlap
merged$region$left <- anno.merged$left
merged$region$right <- anno.merged$right
```


```{R}
all.results <- (data.frame(as.data.frame(merged$region),tabcom, anno.merged,tab.best))
all.results <- all.results[order(all.results$PValue),]
anno.ranges.merged <- detailRanges(merged$region, txdb=TxDb.Mmusculus.UCSC.mm10.knownGene, orgdb=org.Mm.eg.db)
head(anno.ranges.merged$overlap)

#for broad datasets
#for broads
all.results.broad <-(data.frame(as.data.frame(broads),tabbroad, anno.broad))
all.results.broad <-all.results.broad[order(all.results.broad$PValue),]
```
Generate DB .bed file for GO analysis
```{R}
is.sig <- merged$region$FDR <= 0.05
library(rtracklayer)
test <- merged$region[is.sig]
test$score <- -10*log10(merged$region$FDR[is.sig])
names(test) <-paste0("region", 1:sum(is.sig))
export.bed(test, "H3K27ac_4FC_AD12_DBsites_fixedmodel_11072018.bed")
head(read.table("H3K27ac_4FC_AD12_DBsites_fixedmodel_11072018.bed"))
```


Save outputs
```{r}
library(rtracklayer)
save.image(file = "K27ac_12m_4FC_AD_csaw_pseudo_environment_fixed_model_11072018.Rda")
write.csv(all.results, file="H3K27ac_12m_4FC_AD_bcluster_DB_pseudo_all_results_fixedmodel_11072018_test02082018.csv")
write.csv(all.results.broad, file="H3K27ac_12m_4FC_AD_bcluster_DB_pseudo_binned_all_results_broad_fixedmodel_11072018.csv")
```
