for file in *.sam ; do
	
	(samtools view -@62 -b /NGS_Data/Andrew/Alignment/TI_aligned/${file} > /NGS_Data/Andrew/BAM/TI_bam/${file}.bam) 2>/NGS_Data/Andrew/BAM/TI_bam/$file.log

done