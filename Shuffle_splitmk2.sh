for file in *.bam ; do
nlines=$(samtools view ${file} | wc -l) #identify the number of reads in a bam file
nlines=$(( (nlines + 1) / 2 )) #take half of that number
samtools view ${file} | shuf - | split -d -l ${nlines} - "/NGS_Data/Andrew/Pseudoreplicates_BAM/test/${file}_part.sam"
done 