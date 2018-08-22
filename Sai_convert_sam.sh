#Title:sai convert
#Date: 2017-2018
#Author: Andrew Phipps
#notes: calling 2 files for converting sai to sam
#######################################################################################
#Download NGS data from AGRF through wget
for file in *.fq; do
file2="$(basename ${file})_align.sai"
(bwa samse /NGS_Data/Andrew/BWA_alignment/index/mm10.fa $file2 $file > ${file2}.sam ) 2>${file2}.log
done
