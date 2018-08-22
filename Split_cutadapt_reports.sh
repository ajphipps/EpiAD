#!/bin/bash

csplit --digits=2  --quiet --prefix="/NGS_Data/Andrew/Nugen_trimmed_samples/AGRF_K27ac_Nugen_Trimmed/QC/Trimming_reports/" "K27ac_trimming_report.txt" "/This is cutadapt 1.15 with Python 3.5.2/-0" "{*}"

#IFS='This is cutadapt 1.15 with Python 3.5.2' 

#IFS=' '
#read -a TOKEN <<< K27ac_trimming_report.txt
	
#COUNT=1
#for i in "${TOKEN[@]}"; do	# access each element of array
#    	echo IFS > "/NGS_Data/Andrew/test-report/$COUNT.txt"
#    	echo "$i" 
#>> "/NGS_Data/Andrew/test-report/$COUNT.txt"
#   	let "COUNT++"
#done
