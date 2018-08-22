IFS='_' 
for file in *.fq ; do
	read -a TOKEN <<< "$file"
	#echo "${TOKEN[4]}"
 
   cutadapt -j 60 -a "${TOKEN[4]}" -o "./output/$file" "$file"
 
done