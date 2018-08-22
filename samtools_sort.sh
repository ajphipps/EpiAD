#samtools sort to index bam files
for file in *.sam ; do
	(samtools sort -@62 ${file} -o /.${file}_sorted.bam) \
	2>$file_sorted.log #this saves the stderr to directory
done