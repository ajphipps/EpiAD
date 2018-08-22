#deepTools create bigwig - needs to be indexed first!!!!
for file in *.sorted.bam ; do 
bamCoverage --bam ${file} -o /NGS_Data/Andrew/test/${file}.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2652783500 
done