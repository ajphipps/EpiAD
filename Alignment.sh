#!/bin/bash
#Alignment using bowtie2
for file in *.fq ; do
	
	(bowtie2 -p60 --met-file "/NGS_Data/Andrew/Alignment/test/test_output/${file}_metrics"  -x /NGS_Data/Andrew/Alignment/mm10ref/mm10btref -U "$file" -S "/NGS_Data/Andrew/Alignment/K27ac_aligned/$file.sam") 2>/NGS_Data/Andrew/Alignment/K27ac_aligned/$file.log

done