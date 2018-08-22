#script for shuffling and splitting reads in 2
for file in *.bam; do
	echo 'counting number of reads...'
	nlines= $(samtools view -@62 ${file} | wc -l)
	echo nlines
	nlines=$(( (nlines + 1)/2 )) #take half the number
	echo 'shuffling and splitting file in half'
	samtools view $file | shuf - | split -d -l ${nlines} - "NGS_Data/Andrew/Pseudoreplicate_BAM/K27ac/temp_sam/${file}_part.sam"
	done 
