#Author: Andrew Phipps
#Date: 20042018
#notes: summary of DiffBind experiments - K4me3
###############################################################################
#Setting things up:
source("https://bioconductor.org/biocLite.R")
biocLite("DiffBind")
biocLite("edgeR")
biocLite("csaw")
biocLite("DiffBind")
biocLite("RUVSeq")
library(DiffBind)
library(edgeR)
library(rtracklayer)
library(dplyr)
library(RUVSeq)
###############################################################################
#read in sample-sheet for analysis - use unix readlink *.bam to generate filpath table for analysis
samples <- read.csv(file.path("/NGS_Data/Andrew/RStudio/DiffBind/samplesheetK4me3_sorted2.csv"))
#now load in peaks from samples
K4me3_peaks <- dba(sampleSheet=samples)

#plot correlating for binding matrix
plot(K4me3_peaks)
#plot peaks into summits for narrow mark
K4me3_peaks <- dba.count(K4me3_peaks, summits=250)

#now we need to establish a contrast between the samples to test differential binding on 
#eg WT/TG or time-point
K4me3_peaks <- dba.contrast(K4me3_peaks, categories=DBA_CONDITION, block=DBA_TISSUE)

#perform differential analysis
K4me3_peaks <- dba.analyze(K4me3_peaks)

#correlation heatmap can be plotted based on the results of the analysis
plot(K4me3_peaks, contrast=1)

#create a report showing DB sites from K4me3 data
K4me3_DB <- dba.report(K4me3_peaks)

#create a PCA plot for K4me3 analysis
dba.plotPCA(K4me3_peaks,DBA_TISSUE, label=DBA_CONDITION)

#MA plot for showing normalisation  and identifying DB sites
dba.plotMA(K4me3_peaks)

#Volcano plot for showing DB sites - not as useful with 38 total sites
dba.plotVolcano(K4me3_peaks)

corvals <- dba.plotHeatmap(K4me3_peaks, contrast=1, correlations = FALSE)

#compare edgeR and Deseq2 analysis
K4me3_peaks <- dba.analyze(K4me3_peaks, method=DBA_ALL_METHODS)