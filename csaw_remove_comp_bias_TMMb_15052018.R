CSAW removal of COMP bias TMM binned - 15052018
#title: "CSAW_rem_comp_bias_13052018 - 2FC, write output"
#author: "Andrew Phipps"
#date: "13 May 2018"
#output: html_document
#---
#download and install required packages - note you need R 3.3.3 or later
#source("https://bioconductor.org/biocLite.R")
#biocLite("edgeR")
#biocLite("csaw")
#biocLite("DiffBind")
#install.packages("tidyverse")
#install.packages("dplyr")
#biocLite("RUVSeq")
#biocLite("TxDb.Mmusculus.UCSC.mm10.knownGene")
#biocLite("GenomicRanges")

#load required packages into R

library(csaw)
library(edgeR)
library(rtracklayer)
library(dplyr)

library("GenomicRanges")

#import required bam files - for now stick to single timepoint WT vs TG for K4me3 - but not sure if it would be better to load all data
#for multiple marks at a single time point or import all data from the experiment 
#load samples in from marked duplicated files

K4me3_bam.files <- file.path("/NGS_Data/Andrew/Bowtie2_MarkDuplicates/K4me3" , c("512_WT3_K4me3_CC0TAANXX_CATCCAAG_all_trimmed.fq.sam.bam.sorted.bam_markdup.bam", 
                                                                                 "513_WT3_K4me3_CC0TAANXX_GTCAACAG_all_trimmed.fq.sam.bam.sorted.bam_markdup.bam", 
                                                                                 "517_WT3_K4me3_CC0TAANXX_TCGCTATC_all_trimmed.fq.sam.bam.sorted.bam_markdup.bam",
                                                                                 "363_WT3_K4me3_CC0TAANXX_AGCCTATC_all_trimmed.fq.sam.bam.sorted.bam_markdup.bam",
                                                                                 "368_WT3_K4me3_CC0TAANXX_TCGGATTC_all_trimmed.fq.sam.bam.sorted.bam_markdup.bam",
                                                                                 "525_TG3_K4me3_CC0TAANXX_CGGAGTAT_all_trimmed.fq.sam.bam.sorted.bam_markdup.bam",
                                                                                 "534_TG3_K4me3_CC0TAANXX_GAACCTTC_all_trimmed.fq.sam.bam.sorted.bam_markdup.bam",
                                                                                 "535_TG3_K4me3_CC0TAANXX_AGAGGATG_all_trimmed.fq.sam.bam.sorted.bam_markdup.bam",
                                                                                 "374_TG3_K4me3_CC0TAANXX_ACGCTTCT_all_trimmed.fq.sam.bam.sorted.bam_markdup.bam",
                                                                                 "364_TG3_K4me3_CC0TAANXX_CACAGGAA_all_trimmed.fq.sam.bam.sorted.bam_markdup.bam"))
#"157307_TG12_K4me3_CC0TAANXX_ACGTCGTT_all_trimmed.fq.sam.bam.sorted.bam_markdup.bam"))
input <- file.path("/NGS_Data/Andrew/Bowtie2_MarkDuplicates/TI", c("TI_merged.sorted.bam_markdup.bam"))


#design matrix
design <- model.matrix(~factor(c('WT3', 'WT3', 'WT3', 'WT3', 'WT3', 'TG3', 'TG3','TG3','TG3','TG3')))
colnames(design) <- c("intercept", "genotype")

#making blacklist granges with Rtracklayer
library(rtracklayer)
gr_blacklist = import("/NGS_Data/Andrew/bigwig/mm10.blacklist.bed")

###############################################################################
#setting up parameters
frag.length <- 250 #this is the fragment length from sonication
window.width <- 150 #will partition into windows of 150bp - for TF set small window, but for histone marks set window relevent \
##to mark - around 150bp windows - min size of a nucleosome? see 2.5 - Choosing appropriate window size
spacing <- 50
#discard mm10 blacklist regions - remember to include encode blacklist regions #probability mapped incorrectly: 0.01 (mapq 20)
parameters <- readParam(minq=30, discard = gr_blacklist, dedup=TRUE)
data <- windowCounts(K4me3_bam.files, ext=frag.length, width =window.width, param = parameters)

#can visualise the data with 
rowRanges(data)

#visualise counts with 
head(assay(data), n=100)

#RUVseq 
#prelimcontrol <- is.sig.gene
#data1 <- RUVg(data, prelimcontrol, k= 1)

#can visualise to pick fragment size with the following plot - however you need to have marked (not removed) duplicates with Picard
#Sharp spike = fragment length. Low curve = potentially poor IP efficiency (low input from my samples)
max.delay <- 500
dedup.on <- reform(parameters, dedup=TRUE)
plot1 <- correlateReads(K4me3_bam.files, max.delay, param = dedup.on)
plot(0:max.delay, plot1, type='l', ylab="CCF", xlab="Delay (bp)")
#can quantitate the fragment size with the following: 
maximizeCcf(plot1)


#filtering from the 'negative binomial' - log transform of NB is referred to as abundance
##to remove uninteresting regions and lower computational complexity
#a better option for our analysis is to filter via global enrichement
bin.size <- 2000L 
binned1 <- windowCounts(K4me3_bam.files, bin=TRUE, width=bin.size, param=parameters)
filtered.stat <- filterWindows(data, background=binned1, type="global")
#keep samples with a fold change of greater than 3 (3 by default in csaw guide) - try other filtering steps 

filtered.keep <- filtered.stat$filter >log2(3)

#sum of filtered samples
sum(filtered.keep)
#make filtered.data array for downstream analysis
filtered.data <- data[filtered.keep,]

#visualise fold change to confirm that the bulk of background sites are removed by filtering
par(mfrow=c(1,1)) 
hist(filtered.stat$back.abundances, xlab="adjusted bin log-CPM", breaks=100,
     main="", col="grey80", xlim=c(min(filtered.stat$back.abundances), 0))
global.bg <- filtered.stat$abundances -filtered.stat$filter
abline(v=global.bg[1], col="red", lwd=2)
abline(v=global.bg[1]+log2(3), col="blue", lwd=2)
legend("topright", lwd=2, col=c('red', 'blue'), legend=c("Background", "Threshold"))
#Now use global enrichment for elimination of composition bias - here we need to pick one of the two forms of normalisation to be able to correctly quantify

#elimination of composition bias with global enrichment
binned <- windowCounts(K4me3_bam.files, bin=TRUE, width=10000, param=parameters)
filtered.data <- normOffsets(binned, se.out=filtered.data)
filtered.data$norm.factors
#You can test multiple normalisation windows here: too small = low counts and loss of DB, too large and DB will
##be in the same window as background
demo1 <- windowCounts(K4me3_bam.files, bin=TRUE, width=5000, param=parameters)
normOffsets(demo1, se.out=FALSE)
demo2 <- windowCounts(K4me3_bam.files, bin=TRUE, width=15000, param=parameters)
normOffsets(demo, se.out=FALSE)

#visualisation of normalisation with MA plots - generate for each sample 
#vertical shift in the bars might indicate composition bias - ideally want comp factors (line) to pass through
##centre of the cloud
par(mfrow=c(3,3), mar=c(5,4,2,1.5))
adj.counts <- cpm(asDGEList(binned), log=TRUE)
normfacs <- filtered.data$norm.factors
for (i in seq_len(length(K4me3_bam.files)-1)) {
  cur.x <- adj.counts[,1]
  cur.y <- adj.counts[,1+i]
  smoothScatter(x=(cur.x+cur.y)/2+6*log2(10), y=cur.x-cur.y,
                xlab="A", ylab="M", main=paste("1 vs",i+1))
  all.dist <- diff(log2(normfacs[c(i+1, 1)]))
  abline(h=all.dist, col="red")
}

#testing for Diffbinding
#need: filtered.data and filtered.data, original data: K4me3_bam.files, and design matrix
#setting up the data
y<- asDGEList(filtered.data)
#experimental design and setting up data
design
#stabilising estimates with empyrical bayes
y<- estimateDisp(y, design)
summary(y$trended.dispersion)
fit <-glmQLFit(y, design, robust=TRUE)
summary(fit$var.post)


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


#visualise with MDS plots:
#plotMDS()


#Testing for differential binding
results <- glmQLFTest(fit, contrast=c(0,1))
head(results$table)

#assign p-values to co-ordinates
rowData(filtered.data) <-cbind(rowData(filtered.data),results$table)


#examine replicate similarity with MDS plots
par(mfrow=c(2,2), mar=c(5,4,2,2))
adj.counts<-cpm(y,log=TRUE)
for(top in c(100,500,1000,5000)) {
  out <- plotMDS(adj.counts, main=top, col=c("red","red", "red", "red", "red", "blue", "blue", "blue", "blue", "blue" ),
                 labels=c("512", "513", "517","363", "368", "525", "534", "535", "374", "364"), top=top)
}
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
broads <-genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
broads <- resize(broads, width(broads)+3000, fix="end")
head(broads)

#Cluster with external cues/prior information - the mm10 DB

olap <- findOverlaps(broads, rowRanges(filtered.data))
tabbroad <- combineOverlaps(olap, results$table)
head(tabbroad[!is.na(tabbroad$PValue),],n=100)

merged <- mergeWindows(rowRanges(filtered.data), )

is.sig.gene <- tabbroad$FDR <= 0.05
table(tabbroad$direction[is.sig.gene])

tabbroad %>% arrange(PValue) %>% head(n = 10)
tabbroad %>% arrange(FDR) %>% head(n = 100)


#cluster windows into regions without external information/cues - note max.width to limit size (6.2.2)
#tolerance is the min distance for two binding sites to be treated as separate events
mergedwindowsK4me3 <- mergeWindows(rowRanges(filtered.data), tol = 1000, max.width = 2000L)
mergedwindowsK4me3$region
#assigning combined P value for merged windows
p.mergedwindowsK4me3 <- combineTests(mergedwindowsK4me3$id, results$table)
#check to see if most clusters are an acceptable size, if there are huge clusters we need to improve our filtering or limit
summary(width(mergedwindowsK4me3$region))

#now assign direction of fold change to the p value
direction.p.mergedwindowsK4me3 <- p.mergedwindowsK4me3$FDR <= 0.05
table(mergedwindowsK4me3$direction[direction.p.mergedwindowsK4me3])

options(digits = 22)
p.mergedwindowsK4me3 %>% arrange(PValue) %>% head(n = 100)
p.mergedwindowsK4me3 %>% arrange(FDR) %>% head(n = 100)

is.sig.gene.p <- tabbroad$PValue <= 0.05
table(tabbroad$direction[is.sig.gene.p])

#total number of sites that also have acceptable FDR
is.sig.gene.FDR <- tabbroad$FDR <= 0.05
table(tabbroad$direction[is.sig.gene.FDR])

#Write outputs to csv:
write.csv(tabbroad, file= "CSAW_K4me312mtabbroad_mapped_to_genes_markdup_14052018")
write.csv(p.mergedwindowsK4me3, file= "CSAW_K4me312m_basic_window_cluster_markdup_14052018")



