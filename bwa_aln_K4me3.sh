for file in *.fq ; do
(bwa aln -t 62 "/NGS_Data/Andrew/BWA_alignment/index/mm10.fa" ${file} > /NGS_Data/Andrew/BWA_alignment/K4me3_sai/${file}_align.sai) 2> /NGS_Data/Andrew/BWA_alignment/K4me3_sai/${file}.log
done