#Title:R plot venn diagram from bedfiles
#Date: 2017-2018
#Author: Andrew Phipps
#notes: R script to plot venn diagrams from diffBind venn diagrams
#notes:
#######################################################################################
source("https://bioconductor.org/biocLite.R")
#biocLite("edgeR")
#biocLite("csaw")
#biocLite("DiffBind")
#install.packages("tidyverse")
#install.packages("dplyr")
#biocLite("Gviz")
#biocLite("VennDiagrams")
#biocLite("rtracklayer")
Library("rtracklayer")
Library("VennDiagrams")
Library("Gviz")
Library("dplyr")
Library("edgeR")
bed1 <- import.bed(file.path(/NGS_Data/Andrew/Rstudio/csaw)
bed2 <-