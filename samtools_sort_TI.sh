#samtools sort to index bam files
for file in *.bam ; do
	
	(samtools sort -@62 ${file} -o /NGS_Data/Andrew/BAM_sort/TI_bam_sort/${file}.sorted.bam) 2>/NGS_Data/Andrew/BAM_sort/TI_bam_sort/$file.log

done