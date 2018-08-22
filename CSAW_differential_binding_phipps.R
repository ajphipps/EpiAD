#Title: Differential binding analysis with CSAW
#Author: Andrew Phipps
#Date: 06032018
#notes: summary of csaw experiments

###############################################################################
#download and install required packages - note you need R 3.3.3 or later
#source("https://bioconductor.org/biocLite.R")
#biocLite("edgeR")
#biocLite("csaw")
#biocLite("DiffBind")

#load required packages into R
#tutorial link https://www.bioconductor.org/help/course-materials/2015/BioC2015/csaw_lab.html
library(DiffBind)
library(csaw)
library(edgeR)
library(rtracklayer)
library(ggplot2)
#import required bam files - for now stick to single timepoint WT vs TG for K4me3 - but not sure if it would be better to load all data
#for multiple marks at a single time point or import all data from the experiment 
#load samples in from marked duplicated files
K4me3_bam.files <- file.path("/NGS_Data/Andrew/Bowtie2_MarkDuplicates/K4me3" , c("154467_TG12_K4me3_CC0TAANXX_CGTCCATT_all_trimmed.fq.sam.bam.sorted.bam_markdup.bam", 
                                                                         "155542_TG12_K4me3_CC0TAANXX_GTCCTGTT_all_trimmed.fq.sam.bam.sorted.bam_markdup.bam", 
                                                                         "155668_WT12_K4me3_CC0TAANXX_AGCTAGTG_all_trimmed.fq.sam.bam.sorted.bam_markdup.bam",
                                                                         "155669_TG12_K4me3_CC0TAANXX_AGCCGTAA_all_trimmed.fq.sam.bam.sorted.bam_markdup.bam",
                                                                         "155688_WT12_K4me3_CC0TAANXX_CACGTCTA_all_trimmed.fq.sam.bam.sorted.bam_markdup.bam",
                                                                         "155691_TG12_K4me3_CC0TAANXX_GAGTAGAG_all_trimmed.fq.sam.bam.sorted.bam_markdup.bam",
                                                                         "156508_WT12_K4me3_CC0TAANXX_ACTATCGC_all_trimmed.fq.sam.bam.sorted.bam_markdup.bam",
                                                                         "156509_WT12_K4me3_CC0TAANXX_GCGTATCA_all_trimmed.fq.sam.bam.sorted.bam_markdup.bam",
                                                                         "157306_WT12_K4me3_CC0TAANXX_ACTCTCCA_all_trimmed.fq.sam.bam.sorted.bam_markdup.bam",
                                                                         "157307_TG12_K4me3_CC0TAANXX_ACGTCGTT_all_trimmed.fq.sam.bam.sorted.bam_markdup.bam"))

input <- file.path("/NGS_Data/Andrew/Bowtie2_MarkDuplicates/TI", c("TI_merged.sorted.bam_markdup.bam"))
#design matrix
design <- model.matrix(~factor(c('TG12', 'TG12', 'WT12', 'TG12', 'WT12', 'TG12', 'WT12','WT12','WT12','TG12')))
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
parameters <- readParam(minq=30, discard = gr_blacklist)
#parallelisation of 
data <- windowCounts(K4me3_bam.files, ext=frag.length, width =window.width, param = parameters)
#can visualise the data with 
rowRanges(data)
#visualise counts with 
assay(data)
###############################################################################
#can visualise to pick fragment size with the following plot - however you need to have marked (not removed) duplicates with Picard
#Sharp spike = fragment length. Low curve = potentially poor IP efficiency (low input from my samples)
max.delay <- 500
dedup.on <- reform(parameters, dedup=TRUE)
plot1 <- correlateReads(K4me3_bam.files, max.delay, param = dedup.on)
plot(0:max.delay, plot1, type='l', ylab="CCF", xlab="Delay (bp)")
#can quantitate the fragment size with the following: 
maximizeCcf(plot1)
#you can also perform library specific fragmentation - see manual
###############################################################################
###############################################################################
#Filtering steps:
#filtering out low quality reads

#Independent filtering of count data - see Chapter 3
#filtering from the 'negative binomial' - log transform of NB is referred to as abundance
##to remove uninteresting regions and lower computational complexity
#a better option for our analysis is to filter via global enrichement
bin.size <- 2000L 
binned1 <- windowCounts(K4me3_bam.files, bin=TRUE, width=bin.size, param=parameters)
filtered.stat <- filterWindows(data, background=binned1, type="global")
#keep samples with a fold change of greater than 4 (3 by default in csaw guide)
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

###############################################################################
###############################################################################


#elimination of composition bias with global enrichment

binned <- windowCounts(K4me3_bam.files, bin=TRUE, width=10000, param=parameters)
filtered.data <- normOffsets(binned, se.out=filtered.data)
filtered.data$norm.factors
#You can test multiple normalisation windows here: too small = low counts and loss of DB, too large and DB will
##be in the same window as background
#demo <- windowCounts(K4me3_bam.files, bin=TRUE, width=5000, param=parameters)
#normOffsets(demo, se.out=FALSE)
# [1] 0.9748893 1.0295585 0.8987019 1.0386579 1.0815877 0.8709669 0.9466737 1.0718893 1.0981895 1.0167509
#demo <- windowCounts(K4me3_bam.files, bin=TRUE, width=15000, param=parameters)
#normOffsets(demo, se.out=FALSE)
#[1] 0.9847623 1.0302603 0.9183524 1.0549877 1.0909148 0.8883423 0.9719159 1.0686444 1.1166129 0.9051679

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
all.dist <-diff(log2(normfacs[c(i+1, 1)]))
abline(h=all.dist, col="red")
}

###############################################################################
#Eliminating efficiency bias using TMM on high abundance regions 
filtered.data.TMM <- normOffsets(filtered.data, se.out=TRUE)
filtered.data.TMM.efficiency <- filtered.data.TMM$norm.factors
data.comp<- normOffsets(binned, se.out=FALSE)
#visualisation post normalisation 
##Low A-value = background, high A-value = bound. Normalisation factors from removal of comp bias(dashed line), pass 
###through low A-value, removal of efficiency bias pass through (full)
par(mfrow=c(1,2))
bins <- binned
comp <- data.comp
eff <- filtered.data.TMM.efficiency
adjc <-cpm(asDGEList(bins), log=TRUE)
smoothScatter(x=rowMeans(adjc), y=adjc[,1]-adjc[,2], xlab="A-value (background vs whole)", ylab="M", main= "TMM normalisation K4me3 12m")
abline(h=log2(eff[1]/eff[2]), col="red")
abline(h=log2(comp[1]/comp[2]), col="red", lty=2)
###############################################################################
###############################################################################
#testing for Diffbinding
#need: filtered.data.TMM and filtered.data, original data: K4me3_bam.files, and design matrix
#setting up the data
y<- asDGEList(filtered.data.TMM)
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
plotMDS()
#Testing for differential binding
results <- glmQLFTest(fit, contrast=c(0,1))
head(results$table)

#assign p-values to co-ordinates
rowData(filtered.data.TMM) <-cbind(rowData(filtered.data.TMM),results$table)

#examine replicate similarity with MDS plots
par(mfrow=c(2,2), mar=c(5,4,2,2))
adj.counts<-cpm(y,log=TRUE)
for(top in c(100,500,1000,5000)) {
  out <- plotMDS(adj.counts, main=top, col=c("blue","blue", "red", "blue", "red", "blue", "red", "red", "red", "blue" ),
                 labels=c("154467", "155542", "155668","155669", "155688", "155691", "156508", "156509", "157306", "157307"), top=top)
}

###############################################################################
###############################################################################
#Correction for multiple testing
#need - filtered.data.TMM and results
#uses Benamini-Hochbeg method to p-values, which is less conservative than Bonferroni correction, but still provides
##some form of error control. Can correct with regions (wider combinations of windows) - best for broader marks, 
### or use a single window to represent region/cluster - sensible for sharp binding sites

#cluster windows into regions without external information/cues - note max.width to limit size (6.2.2)
#tolerance is the min distance for two binding sites to be treated as separate events
mergedwindowsK4me3 <- mergeWindows(rowRanges(filtered.data.TMM), tol = 100L) #max.width = 8000L)
mergedwindowsK4me3$region
#assigning combined P value for merged windows
p.mergedwindowsK4me3 <- combineTests(mergedwindowsK4me3$id, results$table)
#check to see if most clusters are an acceptable size, if there are huge clusters we need to improve our filtering or limit
summary(width(mergedwindowsK4me3$region))

#now assign direction of fold change to the p value
direction.p.mergedwindowsK4me3 <- p.mergedwindowsK4me3$FDR <= 0.1
table(mergedwindowsK4me3$direction[direction.p.mergedwindowsK4me3])



#option to select only a particular file type from a working directory


##sys.glob("*.bam")
