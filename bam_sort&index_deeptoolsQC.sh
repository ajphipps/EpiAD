#Title: samtools + deeptools for QC
#Author: Andrew Phipps
#Date: 2017-2018
#notes: https://deeptools.readthedocs.io/en/
#notes:

#install samtools 
wget -O /where/to/install https://github.com/samtools/samtools/releases/download/1.7/samtools-1.7.tar.bz2

cd samtools-1.7    # and similarly for bcftools and htslib
./configure --prefix=/where/to/install
make
make install

#add samtools to path
PATH=/where/to/install/bin:$PATH 
export PATH

#samtools sort bam files
for file in *.bam ; do
	(samtools sort -@62 ${file} -o /NGS_Data/Andrew/BAM_sort/${file}.sorted.bam) \
	2>/NGS_Data/Andrew/BAM/$file.log #this saves the stderr to directory
done

#samtools index bam file
for file in *.sorted.bam ; do
	samtools -@62 index $file
done

#note - usually here you would also use picard to markduplicates if you have any for count based data (ChIP-seq or RNA-seq depending on your pipeline)

#deeptools QC for correlation plots or PCA
#install deeptools with python pip
pip install deepTools

#generate summary .npz files from bams
for file in *.bam ; do 
multiBamSummary bins --bamfiles ${file} 
-out ${file}.npz \
--outRawCounts ${file}readCounts_K4me3.tab \
--smartLabels \
--numberOfProcessors 60 \
done

#from here you can generate the PCA or correlation plot (see deeptools website)