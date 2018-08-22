#bamcoverageTest
bamCoverage --bam 525* -o /NGS_Data/Andrew/test/525TItest.bw \
    --binSize 10
    --normalizeUsing RPGC 
    #RPGC (per bin) = number of reads per bin / scaling factor for 1x average coverage. This scaling factor, in turn, is determined from the sequencing depth: (total number of mapped reads * fragment length) / effective genome size. 
    #The scaling factor used is the inverse of the sequencing depth computed for the sample to match the 1x coverage. This option requires –effectiveGenomeSize. 
    #Each read is considered independently, if you want to only count one mate from a pair in paired-end data, then use the –samFlagInclude/–samFlagExclude options.
    --effectiveGenomeSize 2652783500 #normalise to mouse genome size with multimapping - if samples are not multimapped use 2308125349 instead
    --extendReads