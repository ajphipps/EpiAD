#take the 4th string (sequencing tag) and input it into  cutadapt to remove that sequence from the file 
#eg 163734_WT6_K27ac_CC0TAANXX_CATACGGA_all_trimmed.fq
#take CATACGGA and input it as token 4
#!/bin/bash

IFS='_' 
for file in *.fq ; do
	read -a TOKEN <<< "$file"
	#echo "${TOKEN[4]}"
 
   cutadapt -j 60 -a "${TOKEN[4]}" -o "/NGS_Data/Andrew/test/cutadapt/out/$file" "$file"
 
done