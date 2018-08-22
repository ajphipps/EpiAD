for file in *.sam ; do
	
	(samtools view -@62 -b ${file} > /NGS_Data/Andrew/BAM/K4me3_bam/${file}.bam) 2>/NGS_Data/Andrew/BAM/K4me3_bam/$file.log

done