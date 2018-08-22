for file in *dedup.bam ; do
	
	(samtools sort -@62 ${file} -o /NGS_Data/Andrew/Merged_bams/sorted/bam_remove_duplicates/${file}.sorted.bam) 2>/NGS_Data/Andrew/Merged_bams/sorted/bam_remove_duplicates/$file.log

done