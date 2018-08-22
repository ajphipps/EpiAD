#MarkDuplicates for folder
for file in *.sorted.bam ; do	
(java -jar /NGS_Data/Andrew/Tools/picard.jar MarkDuplicates I=${file} O=/NGS_Data/Andrew/Merged_bams/sorted/bam_remove_duplicates/${file}_dedup.bam M=/NGS_Data/Andrew/Merged_bams/sorted/bam_remove_duplicates/${file}_metrics.txt) >2/NGS_Data/Andrew/Merged_bams/sorted/bam_remove_duplicates/
done
