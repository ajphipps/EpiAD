#plain bash way of catting files together
for file in `ls *.fastq | cut -d"_" -f-5`; do
  cat ${file}_* > /NGS_Data/Andrew/AGRF_12978_cat/${file}_all.fastq
done
