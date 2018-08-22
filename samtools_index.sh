#samtools indexing for files in current folder
for file in *.sorted.bam ; do
	samtools index $file
done