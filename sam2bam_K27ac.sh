for file in *.sam ; do
	
	(samtools view -@62 -b /NGS_Data/Andrew/Alignment/K27ac_aligned/${file} > /NGS_Data/Andrew/BAM/K27ac_bam/${file}.bam) 2>/NGS_Data/Andrew/BAM/K27ac_bam/$file.log

done