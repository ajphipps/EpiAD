#Title:BWA_samtools_flagstat
#Date: 2017-2018
#Author: Andrew Phipps
#notes: samtools flagstat to generate alignment  statistics post BWA alignment
#notes:
#######################################################################################
for file in *.bam ; do 
(samtools flagstat -@62 ${file}) 1>./${file}_alnstats.txt
done
 