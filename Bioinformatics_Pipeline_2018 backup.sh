#Title:Bioinformatics_pipeline_2018
#Date: 2017-2018
#Author: Andrew Phipps
#notes: this is a example pathway from NGS data off illumina hiseq 2500 to peak calling and downstream analysis - be sure to read through notes and associated onenote available at: for more information
#notes:

#######################################################################################
#Download NGS data from AGRF through wget

#mount disk for datastorage
sudo mkdir /NGS_Data
sudo mount /dev/vdb /NGS_Data -t auto
#mount tape backup: 
sudo apt-get update
sudo apt-get install nfs-common
#create a mount point on your system**, e.g.
sudo mkdir -p /my/rd/waad
#Add an entry to your File System table (/etc/fstab), use your fav text editor
144.6.253.12:/rd1/waad /my/rd/waad nfs nfsvers=3,rw,nosuid,hard,intr,bg,proto=tcp 0 0
#Back at the command line, ensure your new entry is mounted using
sudo mount --all
#How to check volume space
df -hT
du -h <path/to/directory>
#######################################################################################

#download data
wget -r --no-host-directories -p /NGS_Data/K4me3_K27ac --ftp-user=PipTaberlay --ftp-password=taberlaywoodhouse ftp://melbourne-ftp.agrf.org.au/AGRF_CAGRF13394_CC0TAANXX
#make a copy of the files with cp -a (or R) /folder/where/files/are/saved /folder where you want them
#extract files from new location
gunzip *.gz
#QC data to ID lane effects prior to cat - ensure sample lanes are identical
fastqc .fastq
#cat all files of the same type and sample
#eg each file name is called: 154467_TG12_K27ac_CC0TAANXX_GTGGTATG_L004_R1.fastq
ls 154467******K27ac*.fastq
cat 154467******K27ac*.fastq > 154467_TG12_K27ac_all.fastq
#repeat for each sample 
#ls 154467******K4me3*.fastq
#cat 154467******K4me3*.fastq > 154467K4me3.fastq
#plain bash way of catting files together
for file in `ls *.fastq | cut -d"_" -f-5`; do
  cat ${file}_* > /NGS_Data/Andrew/AGRF_12978_cat/${file}_all.fastq
done

#fastqc each resulting sample and record QC data -t denotes number of threads to run for simultaneous processing
fastqc -t60 *all.fastq -o /NGS_Data/Andrew/CAGRF13394_QC
#also fastqc lane data for future reference and look for lane variation between samples
mkdir /NGS_Data/Andrew/CAGRF13394_QC/laneqc/
#Fastqc v0.11.7
fastqc *R1.fastq -o /NGS_Data/Andrew/CAGRF13394_QC/laneqc/
#store fastqc for further analysis and summary table
#Multiqc the fastqc files to generate summary report for all files version 1.4
Multiqc -o /NGS_Data/Andrew/QualityControl/MultiQC/ -f --interactive ./  

#######################################################################################
#Trimming adapter sequences found from FASTQC - cutadapt or trim_galore - trim_galore is a wrapper for cutadapt
#that automatically trims illumina default adaptors 0.4.3 - updated 07122016 cutadapt 1.15 - it will remove bases with phred score less than 20 and sequences with 
trim_galore -fastqc_args "-t 60" -o /NGS_Data/Andrew/AGRF_13394_trim/ /NGS_Data/Andrew/AGRF_13394_cat/154467_TG12_K27ac_CC0TAANXX_GTGGTATG_all.fastq

for file in *all.fastq ; do
	trim_galore -fastqc_args "-t 60" -o /NGS_Data/Andrew/test/test2/ ${file}
done

#to run these files in parallel
parallel -j 60 trim_galore -fastqc -o /NGS_Data/Andrew/test/parallel_test/ {} ::: *all.fastq

#4 samples from K27me3 did not run, 517,525,534,535 - so re run those files with the following script
for file in *513*.fastq *517*.fastq *525*.fastq *534*.fastq *535*.fastq ;
do
	trim_galore -fastqc_args "-t 60" -o /NGS_Data/Andrew/AGRF_K27me3_Trimmed/ ${file}
done

#trim Nugen adaptors from sequences manually with cutadapt - note that there is no info-file output when running in multicore mode
IFS='_' 
for file in *.fq ; do
	read -a TOKEN <<< "$file"
	#echo "${TOKEN[4]}"
 
   cutadapt -j 60 -a "${TOKEN[4]}" -o "./output/$file" "$file"
 
done
#copy and paste a trimming report for all files from terminal window and then split into separate files
#csplit can split based on a string of text - change prefix to location for files to be generated in
csplit --digits=2  --quiet --prefix="/NGS_Data/Andrew/Nugen_trimmed_samples/AGRF_K27me3_Nugen_Trimmed/K27me3_QC/Trimming_reports/" "K27me3_Nugen_trim_reports.txt" "/This is cutadapt 1.15 with Python 3.5.2/-0" "{*}"
#rename each split file to.txt
for name in /NGS_Data/Andrew/test-report/* ; do mv "$name" "$name"_K27me3_trimming_report.txt ; done
#run fastqc on all files for each set post-trimming
fastqc -t 60 *.fq -o /NGS_Data/Andrew/Nugen_trimmed_samples/AGRF_K27me3_Trimmed/QC
# make multiqc reports for each file and check fastqc result
multiqc -o ./multiqc ./
#######################################################################################
#Alignment scripts - maybe worth trying a couple of different alignment tools to see if we can get a better unique mapped rate - try bowtie1 vs bowtie2

#build bt2 reference from mm10.fa downloaded from UCSC - bowtie2 2.4.3 - updated and added manually to server
bowtie2-build --threads 60 mm10.fa mm10bt2ref
#save these samples in a folder to call on for alignment

#align samples to bowtie2 with mm10 reference genome - output stderr to $file.log and output a metrics file for each sample
#!/bin/bash
#Alignment using bowtie2
for file in *.fq ; do
	
	(bowtie2 -p62 --met-file "/NGS_Data/Andrew/Alignment/K4me3_aligned/${file}_metrics"  -x /NGS_Data/Andrew/Alignment/mm10ref/mm10btref -U "$file" -S "/NGS_Data/Andrew/Alignment/K4me3_aligned/$file.sam") 2>/NGS_Data/Andrew/Alignment/K4me3_aligned/$file.log

done
#12 K4me3 samples wouldn't run through bowtie2 until I rebooted the system, remounted the drive and readded bowtie2 to PATH
PATH=$PATH: /link/to/bowtie2
export PATH
#can pipe bowtie2 into samtools: 
for file in *.fq ; do
	
	(bowtie2 -p62 -x /NGS_Data/Andrew/5Bowtie2Alignment/mm10ref/ -U "$file" | samtools sort -o ${file}.bam -) 2>/NGS_Data/Andrew/RStudio/csaw/csaw_demo/$file.log

done
#can check the log files from each sample to see if the sample completed run with summary support

#######################################################################################
#converting SAM to BAM and peak-calling
#SAMtools view to convert .sam to .bam - samtools version: 1.7 - Do not need to view if we are sorting first
#for file in *.sam ; do
	
#	(samtools view -@62 -b /NGS_Data/Andrew/Alignment/K4me3_aligned/${file} > /NGS_Data/Andrew/BAM/K4me3_bam/${file}.bam) 2>/NGS_Data/Andrew/BAM/K4me3_bam/$file.log

#done
#samtools sort bam files
for file in *.bam ; do
	
	(samtools sort -@62 ${file} -o /NGS_Data/Andrew/BAM_sort/K4me3_bam/${file}.sorted.bam) 2>/NGS_Data/Andrew/BAM/K4me3_bam/$file.log

done
#samtools index bam file
for file in *.sorted.bam ; do
	samtools -@62 index $file
done
#######################################################################################
#For TI samples, use samtools merge to merge then re-sort (if needed) so we have a single merged control for peak calling
ubuntu@nectarbioinformatics-tas:/NGS_Data/Andrew/BAM_sort/TI_bam_sort$ samtools merge -@62 TI_merge2.bam 154467_TG12_TI_CB8LUANXX_GCATAGTC_all_trimmed.fq.sam.bam.sorted.bam 155688_WT12_TI_CB8LUANXX_GCACACAA_all_trimmed.fq.sam.bam.sorted.bam 158027_WT24_TI_CB8LUANXX_TGTGGCTT_all_trimmed.fq.sam.bam.sorted.bam 163734_WT6_TI_CB8LUANXX_AACACGCT_all_trimmed.fq.sam.bam.sorted.bam 326_TG6_TI_CB8LUANXX_TTCACGGA_all_trimmed.fq.sam.bam.sorted.bam 512_WT3_TI_CB8LUANXX_TGCTGTGA_all_trimmed.fq.sam.bam.sorted.bam 525_TG3_TI_CB8LUANXX_CCTCGAAT_all_trimmed.fq.sam.bam.sorted.bam
#run from within the folder

#merging 12m samples for grant prelim data
ubuntu@nectarbioinformatics-tas:/NGS_Data/Andrew/BAM_sort/K4me3_bam_sort$ ls *TG12*.sorted.bam
154467_TG12_K4me3_CC0TAANXX_CGTCCATT_all_trimmed.fq.sam.bam.sorted.bam  155691_TG12_K4me3_CC0TAANXX_GAGTAGAG_all_trimmed.fq.sam.bam.sorted.bam
155542_TG12_K4me3_CC0TAANXX_GTCCTGTT_all_trimmed.fq.sam.bam.sorted.bam  157307_TG12_K4me3_CC0TAANXX_ACGTCGTT_all_trimmed.fq.sam.bam.sorted.bam
155669_TG12_K4me3_CC0TAANXX_AGCCGTAA_all_trimmed.fq.sam.bam.sorted.bam
ubuntu@nectarbioinformatics-tas:/NGS_Data/Andrew/BAM_sort/K4me3_bam_sort$ samtools merge -@62 /NGS_Data/Andrew/Merged_bams/K4me3_TG12_Merge.sorted.bam *TG12*.sorted.bam

#######################################################################################
#At this step - make a copy of the bam files and back up, then also create .bigwig files for visualisation on IGV 

#QC mapping with principal components analysis utilising deepTools multibamsummary > plotPCA


#deepTools analysis and QC - needs to be indexed bam first!!!! - run from within directory - can take a while
#this process can take a long time - doesn't seem to support multithreading
multiBamSummary bins --bamfiles *sorted.bam 
-out readcounts_K4me3.npz --outRawCounts readCounts_K4me3.tab --smartLabels --numberOfProcessors 60

#now do a principal components analysis of the ChIP
#######################################################################################
#deepTools create bigwig - needs to be indexed first!!!! - Example from manual
for file in *.sorted.bam ; do 
bamCoverage --bam ${file} -o /NGS_Data/Andrew/test/${file}.bw \
    --binSize 10
    --normalizeUsing RPGC 
    #RPGC (per bin) = number of reads per bin / scaling factor for 1x average coverage. This scaling factor, in turn, is determined from the sequencing depth: (total number of mapped reads * fragment length) / effective genome size. 
    #The scaling factor used is the inverse of the sequencing depth computed for the sample to match the 1x coverage. This option requires –effectiveGenomeSize. 
    #Each read is considered independently, if you want to only count one mate from a pair in paired-end data, then use the –samFlagInclude/–samFlagExclude options.
    --effectiveGenomeSize 2652783500 #normalise to mouse genome size with multimapping - if samples are not multimapped use 2308125349 instead - length of mappable genome
    --extendReads #needs paired end reads
done
#example without normalisation
#deepTools create bigwig - needs to be indexed first!!!!
for file in *.sorted.bam ; do 
(bamCoverage --bam ${file} -o /NGS_Data/Andrew/bigwig/K4me3_bigwig/${file}.bw --binSize 10 --numberOfProcessors 60) 2>/NGS_Data/Andrew/bigwig/K4me3_bigwig/$file.bw.log
done
#Example extending reads to make smoother plots - each file has been normalised to 1x genome coverage
for file in *.sorted.bam ; do 
(bamCoverage --bam ${file} -o /NGS_Data/Andrew/bigwig/K4me3_bigwig/${file}.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --numberOfProcessors 60 --extendReads 200) 2>/NGS_Data/Andrew/bigwig/K4me3_bigwig/$file.bw.log
done
#example 3 comparing to TI with SES correction from Diaz et al 2012. - this method leads to really messy bw tracks

bamCompare -b1 treatment file -b2 TI file -o /NGS_Data/Andrew/bigwig/testcomparison --scaleFactorsMethod SES 

#######################################################################################
#Peak calling with MACS2
#note: need to figure out if we use the --kep-dup tag, and if we should use --broad for K27me3 = do use broad mode, and broad-cutoff 0.1
#note2 16022018 - Kate suggested we use a range of peak callers and see what similarities come out from trying a few options for peak callers - like peak ranger, SICER, MACS2
for file in *.sorted.bam ; do
(macs2 callpeak -g 1.87e9 -q 0.05 -t ${file} -c <path/to/TIs or controls> -n ${file} --outdir /NGS_Data/Andrew/Peak_calling/test ) 2> <path/to/output>
done
#example parallelisation of MACS2
### call broad peaks with models
#got it to work - this one will remove PCR duplicates prior to peak calling, then you can compare fold change between samples and Q value
cat *.sorted.bam | parallel --max-procs=12 'macs2 callpeak -t {} \
 -c {}-G-NC.sorted.bam --broad -g hs --broad-cutoff 0.1 -n {}-A-NC-broad-model -q 0.01 \
 --outdir {}-A-NC-broad-model-peaks 2> {}-A-NC-broad-model.stderr' 
#tests for merged data
#broad for K27me3 - tmpfs issues: https://www.biostars.org/p/152765/
ubuntu@nectarbioinformatics-tas:/NGS_Data/Andrew/Merged_bams/sorted$ cat K27me3samplenames.txt | TMPDIR=/NGS_Data/tmp parallel -j 2 "(macs2 callpeak -g 1.87e9 -q 0.05 -t {} -c ./TI/TI_merged.sorted.bam.sorted.bam --broad -n {}_broad --outdir ./test_peakcall/ ) 2>/NGS_Data/Andrew/Merged_bams/sorted/test_peakcall/{}_broad.log"
#broad relaxed K27me3 - q value the same as Gjoneska (broad mode, relaxed p=0.1)
cat K27me3samplenames.txt | TMPDIR=/NGS_Data/tmp parallel -j 2 "(macs2 callpeak -g 1.87e9 -q 0.01 -t {} -c ./TI/TI_merged.sorted.bam.sorted.bam --broad -n {}_broadrelax --outdir ./test_peakcall/ ) 2>/NGS_Data/Andrew/Merged_bams/sorted/test_peakcall/{}_broadrelax.log"

#also try K27me3 with p= 0.01 and broad-cutoff - compare peaks for each with Diffbind in R
cat K27me3samplenames.txt | TMPDIR=/NGS_Data/tmp parallel -j 2 "(macs2 callpeak -g 1.87e9 -p 0.01 -t {} -c ./TI/TI_merged.sorted.bam.sorted.bam --broad --broad-cutoff 0.1 -n {}_broadrelax_P --outdir ./test_peakcall/ ) 2>/NGS_Data/Andrew/Merged_bams/sorted/test_peakcall/{}_broadrelax_p.log"

#######################################################################################
#K4me3 MACS2 peak calling - call within merged, sorted bam folder
cat K4me3samples.txt | TMPDIR=/NGS_Data/tmp parallel -j 35 "(macs2 callpeak -g 1.87e9 -q 0.05 -t {} -c /NGS_Data/Andrew/BAM_sort/TI_bam_sort/*.sorted.bam -n {} --outdir /NGS_Data/Andrew/Peak_calling/MACS_K4me3/ ) 2>/NGS_Data/Andrew/Peak_calling/MACS_K4me3/{}.log"
#K27ac MACS2 peak calling
cat K27acsamples.txt | TMPDIR=/NGS_Data/tmp parallel -j 35 "(macs2 callpeak -g 1.87e9 -q 0.05 -t {} -c /NGS_Data/Andrew/BAM_sort/TI_bam_sort/*.sorted.bam -n {} --outdir /NGS_Data/Andrew/Peak_calling/MACS_K27ac/ ) 2>/NGS_Data/Andrew/Peak_calling/MACS_K27ac/{}.log"
#TI MACS2 peak calling against eachother - determine if there is anything we need to exclude from further analysis
cat TIsamples.txt | TMPDIR=/NGS_Data/tmp parallel -j 7 "(macs2 callpeak -g 1.87e9 -q 0.05 -t {} -c /NGS_Data/Andrew/BAM_sort/TI_bam_sort/Merged_bam/TI_merged.sorted.bam -n {} --outdir /NGS_Data/Andrew/Peak_calling/MACS_TI_comparison/ ) 2>/NGS_Data/Andrew/Peak_calling/MACS_TI_comparison/{}.log"
#the output of the TI cat was a LARGE range of peaks 
#MACS2 peakcalling K27me3 with relaxed q value, and broadcutoff 0.1
cat K27me3samplenames.txt | TMPDIR=/NGS_Data/tmp parallel -j 35 "(macs2 callpeak -g 1.87e9 -q 0.01 --broad --broad-cutoff 0.1  -t {} -c /NGS_Data/Andrew/BAM_sort/TI_bam_sort/*.sorted.bam -n {} --outdir /NGS_Data/Andrew/Peak_calling/MACS_K27me3_broadq01/ ) 2>/NGS_Data/Andrew/Peak_calling/MACS_K27me3_broadq01/{}.log"
#tried peak calling with Gjoneska's p=0.01 and it had a over 40K peaks, lots of which had p=0 - do not use their relaxed P value, maybe relaxed Q will work.
cat K27me3samplenames.txt | TMPIR=/NGS_Data/tmp parallel -j 2 "(mac2 callpeak -g 1.87e9 -p 0.01 -t K27me3_TG12_Merge.bam.sorted.bam -c ./TI/TI_merged.sorted.bam.sorted.bam --broad --broad-cutoff 0.1 -n K27me3_TG12_Merge.bam.sorted.bam_broadrelax_P --outdir ./test_peakcall/)"
#current paper - https://www.sciencedirect.com/science/article/pii/S0969996118300391#bb0205 used q 0.1 for cell line work K27ac - test to see what that looks like
cat samplenames.txt | TMPDIR=/NGS_Data/tmp parallel -j 35 "(macs2 callpeak -g 1.87e9 -q 0.01 -t {} -c /NGS_Data/Andrew/Bowtie2BAM_sort/TI_bam_sort/*.sorted.bam -n {}_q01 --outdir /NGS_Data/Andrew/Merged_bams/sorted/test_peakcall/ ) 2>/NGS_Data/Andrew/Merged_bams/sorted/test_peakcall/{}.log"

#######################################################################################
#22022018 - trying alignment with BWA to see if there is much difference between that and bowtie2 - test 2 samples only
#bwa test
#build bwa index from within folder containing mm10.fa - if moving resulting files, include mm10.fa
bwa index mm10.fa
#generates *.64.* alignment files that are used for downstream
#align samples to reference generated with index
bwa aln /NGS_Data/Andrew/BWA_alignment/index/mm10.fa <sample> | 
#convert sai to sam
for file in *.fq; do
file2="$(basename ${file})_align.sai"
(bwa samse /NGS_Data/Andrew/BWA_alignment/index/mm10.fa $file2 $file > ${file2}.sam ) 2>${file2}.log
done
#convert sam to bam
samtools view - | samtools sort - -o /NGS_Data/BWA_alignment/K4me3_sorted_bam
#convert sam to sorted bam 
#samtools sort bam files
for file in *.bam ; do
	(samtools sort -@62 ${file} -o /.${file}_sorted.bam) \
	2>$file_sorted.log #this saves the stderr to directory
done
#index bam

#######################################################################################
#test alignment with Bowtie1
bowtie-build --threads 60 -f mm10.fa bowtie1_mm10

#######################################################################################
#install picard 2.17.10, and MarkDuplicates for all files prior to importing into CSAW - this took a whole day to process all samples. Not possible to run in parallel
## due to amount of memory required per sample - just set loop per ChIP mark
#mark duplicates from merged bams and output to a new location
for file in *.bam ; do	
(java -jar /NGS_Data/Andrew/Tools/picard.jar MarkDuplicates \
	I=${file} \
	O=./${file}_markdup.bam \
	M=./${file}_metrics.txt \
	VALIDATION_STRINGENCY=LENIENT )
2>./${file}_mardup.log
done

cat K4me3samples.txt | TMPDIR=/NGS_Data/tmp parallel -j 35 "( java -jar /NGS_Data/Andrew/Tools/picard.jar MarkDuplicates \
I={} \
O=/NGS_Data/Andrew/Bowtie2_MarkDuplicates/K4me3/{}_markdup.bam \
M=/NGS_Data/Andrew/Bowtie2_MarkDuplicates/K4me3/{}_metrics.txt \
VALIDATION_STRINGENCY=LENIENT ) \
2>/NGS_Data/Andrew/Bowtie2_MarkDuplicates/K4me3/{}_markdup.log "
#######################################################################################
#csaw making reference peakset for clustering
samtools merge -@62 /NGS_Data/Andrew/Bowtie2_MarkDuplicates/K4me3/merged12m/merged_12m_samples.sorted.bam *12m*.bam 
samtools sort -@62  merged_12m_samples.sorted.bam 

#######################################################################################
#DEEPTOOLS ANALYSIS NOTES
#23042018 - looking back at QC samples - note samples that had above 80% duplication from K4me3 and TIs and remove to see if that improves analysis
#157307 K4me3
#163734 WT6 TI 
#readlink *.bam will list full file-path for all samples to easily make the samplesheet for diffbind
#noticed that when calling peaks against the new merged control it hasn't performed correctly. Need to recall peaks with MACS2, while removing the dodgy TI from the sample pool
cat K4me3samples.txt | TMPDIR=/NGS_Data/tmp parallel -j 35 "(macs2 callpeak -g 1.87e9 -q 0.05 -t {} -c /NGS_Data/Andrew/BAM_sort/TI_bam_sort/*.sorted.bam -n {} --outdir /NGS_Data/Andrew/Peak_calling/MACS_K4me3_TI-154467/ ) \
2>/NGS_Data/Andrew/Peak_calling/MACS_K4me3_TI-154467/{}.log"

#######################################################################################
#14052018 
#re-peak call K4me3, K27me3 and K27ac samples without the noisy 163734 TI input control
#Also going back to -q=0.05 to stay similar to K4me3 and K27ac samples, but still include the --broad and --broad-cutoff 0.1 modifier
cat K27me3samplenames.txt | TMPDIR=/NGS_Data/tmp parallel -j 35 \
"(macs2 callpeak -g 1.87e9 -q 0.05 --broad --broad-cutoff 0.1 -t {} -c /NGS_Data/Andrew/Bowtie2BAM_sort/TI_bam_sort/*.sorted.bam -n {} --outdir /NGS_Data/Andrew/Peak_calling/MACS_K27me3_TI-163734/ ) \2>/NGS_Data/Andrew/Peak_calling/MACS_K27me3_TI-163734/{}.log"



#######################################################################################
#date:08062018
#merge biological replicates and make pseudo-replicates for my samples
samtools merge -u -@62 /NGS_Data/Andrew/Pseudoreplicate_BAM/H3K27me3_intermediate/H3K27me3_WT3_merged.bam *WT3*.bam 

#script to shuffle and split merged bam files
for file in *.bam ; do
echo ""$file" "is being processed...""
samtools view -H ${file} > ${file}_header.sam 
nlines=$(samtools view -@62 ${file} | wc -l) #identify the number of reads in a bam file
echo "total number of reads:"
echo "$nlines" 
nlines=$(( (nlines + 1) / 2 )) #take half of that number
samtools view -@62 ${file} | shuf - | split -d -l ${nlines} - "/NGS_Data/Andrew/Pseudoreplicate_BAM/H3K27ac_intermediate/test/${file}_pseudo.sam"
done 

#this will give a header file, and 2 pseudoreplicates that need to be catted together - at this stage manually
cat K27ac_WT3_merged.bam_header.sam K27ac_WT3_merged.bam_pseudo.sam00 | samtools view -bS -@62 - | samtools sort -@62 - -o ../H3K27ac_pseudoreplicates/K27ac_WT3_pseudo1.sorted.bam &
cat K27ac_WT3_merged.bam_header.sam K27ac_WT3_merged.bam_pseudo.sam01 | samtools view -bS -@62 - | samtools sort -@62 - -o ../H3K27ac_pseudoreplicates/K27ac_WT3_pseudo2.sorted.bam &
cat K27ac_TG3_merged.bam_header.sam K27ac_TG3_merged.bam_pseudo.sam00 | samtools view -bS -@62 - | samtools sort -@62 - -o ../H3K27ac_pseudoreplicates/K27ac_TG3_pseudo1.sorted.bam &
cat K27ac_TG3_merged.bam_header.sam K27ac_TG3_merged.bam_pseudo.sam01 | samtools view -bS -@62 - | samtools sort -@62 - -o ../H3K27ac_pseudoreplicates/K27ac_TG3_pseudo2.sorted.bam &
cat K27ac_TG6_merged.bam_header.sam K27ac_TG6_merged.bam_pseudo.sam01 | samtools view -bS -@62 - | samtools sort -@62 - -o ../H3K27ac_pseudoreplicates/K27ac_TG6_pseudo2.sorted.bam &
cat K27ac_WT6_merged.bam_header.sam K27ac_WT6_merged.bam_pseudo.sam00 | samtools view -bS -@62 - | samtools sort -@62 - -o ../H3K27ac_pseudoreplicates/K27ac_WT6_pseudo1.sorted.bam &
cat K27ac_WT6_merged.bam_header.sam K27ac_WT6_merged.bam_pseudo.sam01 | samtools view -bS -@62 - | samtools sort -@62 - -o ../H3K27ac_pseudoreplicates/K27ac_WT6_pseudo2.sorted.bam &
cat K27ac_TG6_merged.bam_header.sam K27ac_TG6_merged.bam_pseudo.sam00 | samtools view -bS -@62 - | samtools sort -@62 - -o ../H3K27ac_pseudoreplicates/K27ac_TG6_pseudo1.sorted.bam &
cat K27ac_TG12_merged.bam_header.sam K27ac_TG12_merged.bam_pseudo.sam01 | samtools view -bS -@62 - | samtools sort -@62 - -o ../H3K27ac_pseudoreplicates/K27ac_TG12_pseudo2.sorted.bam &
cat K27ac_WT12_merged.bam_header.sam K27ac_WT12_merged.bam_pseudo.sam00 | samtools view -bS -@62 - | samtools sort -@62 - -o ../H3K27ac_pseudoreplicates/K27ac_WT12_pseudo1.sorted.bam &
cat K27ac_WT12_merged.bam_header.sam K27ac_WT12_merged.bam_pseudo.sam01 | samtools view -bS -@62 - | samtools sort -@62 - -o ../H3K27ac_pseudoreplicates/K27ac_WT12_pseudo2.sorted.bam &
cat K27ac_TG12_merged.bam_header.sam K27ac_TG12_merged.bam_pseudo.sam00 | samtools view -bS -@62 - | samtools sort -@62 - -o ../H3K27ac_pseudoreplicates/K27ac_TG12_pseudo1.sorted.bam &
cat K27ac_WT24_merged.bam_header.sam K27ac_WT24_merged.bam_pseudo.sam00 | samtools view -bS -@62 - | samtools sort -@62 - -o ../H3K27ac_pseudoreplicates/K27ac_WT24_pseudo1.sorted.bam &
cat K27ac_WT24_merged.bam_header.sam K27ac_WT24_merged.bam_pseudo.sam01 | samtools view -bS -@62 - | samtools sort -@62 - -o ../H3K27ac_pseudoreplicates/K27ac_WT24_pseudo2.sorted.bam &
-----------
cat H3K27me3_TG12_merged.bam_header.sam H3K27me3_TG12_merged.bam_pseudo.sam00 | samtools view -bS -@62 - | samtools sort -@62 - -o ../H3K27me3_pseudoreplicates/K27me3_TG12_pseudo1.sorted.bam &
cat H3K27me3_TG12_merged.bam_header.sam H3K27me3_TG12_merged.bam_pseudo.sam01 | samtools view -bS -@62 - | samtools sort -@62 - -o ../H3K27me3_pseudoreplicates/K27me3_TG12_pseudo2.sorted.bam &
cat H3K27me3_WT12_merged.bam_header.sam H3K27me3_WT12_merged.bam_pseudo.sam00 | samtools view -bS -@62 - | samtools sort -@62 - -o ../H3K27me3_pseudoreplicates/K27me3_WT12_pseudo1.sorted.bam &
cat H3K27me3_WT12_merged.bam_header.sam H3K27me3_WT12_merged.bam_pseudo.sam01 | samtools view -bS -@62 - | samtools sort -@62 - -o ../H3K27me3_pseudoreplicates/K27me3_WT12_pseudo2.sorted.bam &
cat H3K27me3_TG6_merged.bam_header.sam H3K27me3_TG6_merged.bam_pseudo.sam00 | samtools view -bS -@62 - | samtools sort -@62 - -o ../H3K27me3_pseudoreplicates/K27me3_TG6_pseudo1.sorted.bam &
cat H3K27me3_TG6_merged.bam_header.sam H3K27me3_TG6_merged.bam_pseudo.sam01 | samtools view -bS -@62 - | samtools sort -@62 - -o ../H3K27me3_pseudoreplicates/K27me3_TG6_pseudo2.sorted.bam &
cat H3K27me3_WT6_merged.bam_header.sam H3K27me3_WT6_merged.bam_pseudo.sam00 | samtools view -bS -@62 - | samtools sort -@62 - -o ../H3K27me3_pseudoreplicates/K27me3_WT6_pseudo1.sorted.bam &
cat H3K27me3_WT6_merged.bam_header.sam H3K27me3_WT6_merged.bam_pseudo.sam01 | samtools view -bS -@62 - | samtools sort -@62 - -o ../H3K27me3_pseudoreplicates/K27me3_WT6_pseudo2.sorted.bam &
cat H3K27me3_TG3_merged.bam_header.sam H3K27me3_TG3_merged.bam_pseudo.sam00 | samtools view -bS -@62 - | samtools sort -@62 - -o ../H3K27me3_pseudoreplicates/K27me3_TG3_pseudo1.sorted.bam &
cat H3K27me3_TG3_merged.bam_header.sam H3K27me3_TG3_merged.bam_pseudo.sam01 | samtools view -bS -@62 - | samtools sort -@62 - -o ../H3K27me3_pseudoreplicates/K27me3_TG3_pseudo2.sorted.bam &
cat H3K27me3_WT3_merged.bam_header.sam H3K27me3_WT3_merged.bam_pseudo.sam00 | samtools view -bS -@62 - | samtools sort -@62 - -o ../H3K27me3_pseudoreplicates/K27me3_WT3_pseudo1.sorted.bam &
cat H3K27me3_WT3_merged.bam_header.sam H3K27me3_WT3_merged.bam_pseudo.sam01 | samtools view -bS -@62 - | samtools sort -@62 - -o ../H3K27me3_pseudoreplicates/K27me3_WT3_pseudo2.sorted.bam &
cat H3K27me3_WT24_merged.bam_header.sam H3K27me3_WT24_merged.bam_pseudo.sam00 | samtools view -bS -@62 - | samtools sort -@62 - -o ../H3K27me3_pseudoreplicates/K27me3_WT24_pseudo1.sorted.bam &
cat H3K27me3_WT24_merged.bam_header.sam H3K27me3_WT24_merged.bam_pseudo.sam01 | samtools view -bS -@62 - | samtools sort -@62 - -o ../H3K27me3_pseudoreplicates/K27me3_WT24_pseudo2.sorted.bam &
-----------
cat H3K4me3_TG12_merged.bam_header.sam H3K4me3_TG12_merged.bam_pseudo.sam00 | samtools view -bS -@62 - | samtools sort -@62 - -o ../H3K4me3_pseudoreplicates/H3K4me3_TG12_pseudo1.sorted.bam &
cat H3K4me3_TG12_merged.bam_header.sam H3K4me3_TG12_merged.bam_pseudo.sam01 | samtools view -bS -@62 - | samtools sort -@62 - -o ../H3K4me3_pseudoreplicates/H3K4me3_TG12_pseudo2.sorted.bam &
cat H3K4me3_WT12_merged.bam_header.sam H3K4me3_WT12_merged.bam_pseudo.sam00 | samtools view -bS -@62 - | samtools sort -@62 - -o ../H3K4me3_pseudoreplicates/H3K4me3_WT12_pseudo1.sorted.bam &
cat H3K4me3_WT12_merged.bam_header.sam H3K4me3_WT12_merged.bam_pseudo.sam01 | samtools view -bS -@62 - | samtools sort -@62 - -o ../H3K4me3_pseudoreplicates/H3K4me3_WT12_pseudo2.sorted.bam &
cat H3K4me3_TG6_merged.bam_header.sam H3K4me3_TG6_merged.bam_pseudo.sam00 | samtools view -bS -@62 - | samtools sort -@62 - -o ../H3K4me3_pseudoreplicates/H3K4me3_TG6_pseudo1.sorted.bam &
cat H3K4me3_TG6_merged.bam_header.sam H3K4me3_TG6_merged.bam_pseudo.sam01 | samtools view -bS -@62 - | samtools sort -@62 - -o ../H3K4me3_pseudoreplicates/H3K4me3_TG6_pseudo2.sorted.bam &
cat H3K4me3_WT6_merged.bam_header.sam H3K4me3_WT6_merged.bam_pseudo.sam00 | samtools view -bS -@62 - | samtools sort -@62 - -o ../H3K4me3_pseudoreplicates/H3K4me3_WT6_pseudo1.sorted.bam &
cat H3K4me3_WT6_merged.bam_header.sam H3K4me3_WT6_merged.bam_pseudo.sam01 | samtools view -bS -@62 - | samtools sort -@62 - -o ../H3K4me3_pseudoreplicates/H3K4me3_WT6_pseudo2.sorted.bam &
cat H3K4me3_TG3_merged.bam_header.sam H3K4me3_TG3_merged.bam_pseudo.sam00 | samtools view -bS -@62 - | samtools sort -@62 - -o ../H3K4me3_pseudoreplicates/H3K4me3_TG3_pseudo1.sorted.bam &
cat H3K4me3_TG3_merged.bam_header.sam H3K4me3_TG3_merged.bam_pseudo.sam01 | samtools view -bS -@62 - | samtools sort -@62 - -o ../H3K4me3_pseudoreplicates/H3K4me3_TG3_pseudo2.sorted.bam &
cat H3K4me3_WT3_merged.bam_header.sam H3K4me3_WT3_merged.bam_pseudo.sam00 | samtools view -bS -@62 - | samtools sort -@62 - -o ../H3K4me3_pseudoreplicates/H3K4me3_WT3_pseudo1.sorted.bam &
cat H3K4me3_WT3_merged.bam_header.sam H3K4me3_WT3_merged.bam_pseudo.sam01 | samtools view -bS -@62 - | samtools sort -@62 - -o ../H3K4me3_pseudoreplicates/H3K4me3_WT3_pseudo2.sorted.bam &
cat H3K4me3_WT24_merged.bam_header.sam H3K4me3_WT24_merged.bam_pseudo.sam00 | samtools view -bS -@62 - | samtools sort -@62 - -o ../H3K4me3_pseudoreplicates/H3K4me3_WT24_pseudo1.sorted.bam &
cat H3K4me3_WT24_merged.bam_header.sam H3K4me3_WT24_merged.bam_pseudo.sam01 | samtools view -bS -@62 - | samtools sort -@62 - -o ../H3K4me3_pseudoreplicates/H3K4me3_WT24_pseudo2.sorted.bam &
-----------

#now this file can be indexed
samtools index

#Then perform peak calling as previously - Remember not to include 163734 TI as it failed QC
cat K27ac_samples.txt | TMPDIR=/NGS_Data/tmp parallel -j 35 "(macs2 callpeak -g 1.87e9 -q 0.05 -t {} -c /NGS_Data/Andrew/Bowtie2BAM_sort/TI_bam_sort/*.sorted.bam -n {} --outdir /NGS_Data/Andrew/Pseudoreplicate_BAM/H3K27ac_pseudo_peaks/ ) \
2>/NGS_Data/Andrew/Pseudoreplicate_BAM/H3K27ac_pseudo_peaks/{}.log"
#H3K27me3
cat H3K27me3_samples.txt | TMPDIR=/NGS_Data/tmp parallel -j 35 "(macs2 callpeak -g 1.87e9 -q 0.05 --broad --broad-cutoff 0.1 -t {} -c /NGS_Data/Andrew/Bowtie2BAM_sort/TI_bam_sort/*.sorted.bam -n {} --outdir /NGS_Data/Andrew/Pseudoreplicate_BAM/H3K27me3_pseudoreplicate_peaks/ ) \
2>/NGS_Data/Andrew/Pseudoreplicate_BAM/H3K27me3_pseudoreplicate_peaks/{}.log"
#H3K4me3
cat H3K4me3_samples.txt | TMPDIR=/NGS_Data/tmp parallel -j 35 "(macs2 callpeak -g 1.87e9 -q 0.05 -t {} -c /NGS_Data/Andrew/Bowtie2BAM_sort/TI_bam_sort/*.sorted.bam -n {} --outdir /NGS_Data/Andrew/Pseudoreplicate_BAM/H3K4me3_pseudoreplicate_peaks/ ) \
2>/NGS_Data/Andrew/Pseudoreplicate_BAM/H3K4me3_pseudoreplicate_peaks/{}.log"

#now we can create bigwigs from the pseudoreplicate bam files 
for file in *.sorted.bam ; do 
(bamCoverage --bam ${file} -o /NGS_Data/Andrew/Pseudoreplicate_BAM/bigwig/H3K27ac_bw/${file}.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --numberOfProcessors 60 --extendReads 200) 2>/NGS_Data/Andrew/Pseudoreplicate_BAM/bigwig/H3K27ac_bw/${file}.bw.log
done

for file in *.sorted.bam ; do 
(bamCoverage --bam ${file} -o /NGS_Data/Andrew/Pseudoreplicate_BAM/bigwig/H3K27me3_bw/${file}.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --numberOfProcessors 60 --extendReads 200) 2>/NGS_Data/Andrew/Pseudoreplicate_BAM/bigwig/H3K27me3_bw/${file}.bw.log
done

for file in *.sorted.bam ; do 
(bamCoverage --bam ${file} -o /NGS_Data/Andrew/Pseudoreplicate_BAM/bigwig/H3K4me3_bw/${file}.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --numberOfProcessors 60 --extendReads 200) 2>/NGS_Data/Andrew/Pseudoreplicate_BAM/bigwig/H3K4me3_bw/${file}.bw.log
done
#######################################################################################
#THOR H3K27ac ANALYSIS ATTEMPT WITH ALL BIOLOGICAL REPLICATES
#06072018
#CONFIG FILE
#rep1
/NGS_Data/Andrew/Bowtie2BAM_sort/K27ac_bam_sort/155668_WT12_K27ac_CC0TAANXX_TCACTCGA_all_trimmed.fq.sam.bam.sorted.bam
/NGS_Data/Andrew/Bowtie2BAM_sort/K27ac_bam_sort/155688_WT12_K27ac_CC0TAANXX_CTGTGGTA_all_trimmed.fq.sam.bam.sorted.bam
/NGS_Data/Andrew/Bowtie2BAM_sort/K27ac_bam_sort/156508_WT12_K27ac_CC0TAANXX_ACTCCTAC_all_trimmed.fq.sam.bam.sorted.bam
/NGS_Data/Andrew/Bowtie2BAM_sort/K27ac_bam_sort/156509_WT12_K27ac_CC0TAANXX_CCACAACA_all_trimmed.fq.sam.bam.sorted.bam
/NGS_Data/Andrew/Bowtie2BAM_sort/K27ac_bam_sort/157306_WT12_K27ac_CC0TAANXX_CCGCTTAA_all_trimmed.fq.sam.bam.sorted.bam
#rep2
/NGS_Data/Andrew/Bowtie2BAM_sort/K27ac_bam_sort/154467_TG12_K27ac_CC0TAANXX_GTGGTATG_all_trimmed.fq.sam.bam.sorted.bam
/NGS_Data/Andrew/Bowtie2BAM_sort/K27ac_bam_sort/155542_TG12_K27ac_CC0TAANXX_AACACCAC_all_trimmed.fq.sam.bam.sorted.bam
/NGS_Data/Andrew/Bowtie2BAM_sort/K27ac_bam_sort/155669_TG12_K27ac_CC0TAANXX_GGTGTACA_all_trimmed.fq.sam.bam.sorted.bam
/NGS_Data/Andrew/Bowtie2BAM_sort/K27ac_bam_sort/155691_TG12_K27ac_CC0TAANXX_TCTAGGAG_all_trimmed.fq.sam.bam.sorted.bam
/NGS_Data/Andrew/Bowtie2BAM_sort/K27ac_bam_sort/157307_TG12_K27ac_CC0TAANXX_TGGAAGCA_all_trimmed.fq.sam.bam.sorted.bam
#genome
/NGS_Data/Andrew/mm10/mm10.fa
#chrom_sizes
/NGS_Data/Andrew/THOR/mm10.chrom.sizes
#inputs1
/NGS_Data/Andrew/Bowtie2BAM_sort/TI_bam_sort/Merged_bam/TI_merge_without163734.sorted.bam
#inputs2
/NGS_Data/Andrew/Bowtie2BAM_sort/TI_bam_sort/Merged_bam/TI_merge_without163734.sorted.bam
#######################################################################################
#analysis script
rgt-THOR H3K27ac_AD12m_THOR.config --name H3K27ac_AD12m_test --merge 1000 --pvalue 0.05 --binsize 150 --step 50 --exts 250 --report --output-dir /NGS_Data/Andrew/THOR/H3K27ac_THOR_all --deadzones /NGS_Data/Andrew/THOR/mm10.blacklist.bed  


#######################################################################################
#09072018 - markdup H3K27me3 pseudoreplicates for future analysis
for file in *WT*.bam ; do	
(java -jar /NGS_Data/Andrew/Tools/picard.jar MarkDuplicates \
	I=${file} \
	O=./${file}_markdup.bam \
	M=./${file}_metrics.txt \
	VALIDATION_STRINGENCY=LENIENT )
2>./${file}_mardup.log
done
#markdup didn't work with the K27me3 files - removed the vast majority of sites across the genome and resulted in issues with analysis
#######################################################################################
#DeepTools Heatmap generation
computeMatrix reference-point -S /NGS_Data/Andrew/Pseudoreplicate_BAM/bigwig/H3K27ac_bw/*WT12*.bw /NGS_Data/Andrew/Pseudoreplicate_BAM/bigwig/H3K27ac_bw/*WT24*.bw -R /NGS_Data/Andrew/mm10_knowngene.bed -o ./K27ac_12v24M_TSS_pseudo_matrix -a 5000 -b 5000 -p max -bl /NGS_Data/Andrew/bigwig/mm10.blacklist.bed

#have to make sure the bed file doesn't have any headers etc otherwise computeMatrix will not work
for file in *matrix ; do 
plotHeatmap --matrixFile ${file} --zMin 0 --zMax 7.5  --outFileName ${file}_heatmap.eps --dpi 600 --plotFileFormat eps
done
#eg
computeMatrix reference-point -S /NGS_Data/Andrew/Pseudoreplicate_BAM/bigwig/H3K27ac_bw/*WT3*.bw /NGS_Data/Andrew/Pseudoreplicate_BAM/bigwig/H3K27ac_bw/*WT6*.bw -R /NGS_Data/Andrew/mm10_knowngene.bed -o ./K27ac_3v6M_TSS_pseudo_matrix -a 5000 -b 5000 -p max -bl /NGS_Data/Andrew/bigwig/mm10.blacklist.bed

#can also use the fantom5 enhancer track to compare H3K27ac bigwigs with known mm10 enhancers 
#download fantom5 mm10 enhancer bed file then convert to bed6 format

#######################################################################################
#breadth of coverage test for dead depth 
#http://www.metagenomics.wiki/tools/samtools/breadth-of-coverage
PATH=$PATH:/NGS_Data/Andrew/Tools/bowtie2-2.3.4-linux-x86_64/
export PATH

#small loop for calculating the % genome covered by reads 
for file in *.bam
samtools mpileup ${file} | awk -v X="${5}" '$4>=5' | wc -l > K27ac_breadth.txt
done
#total number of reads in reference genome
bowtie2-inspect-s /NGS_Data/Andrew/mm10/mm10btref | awk '-F "\t" BEGIN{L=0}; {L=L+$3}; END{print L}'
#after this you can calculate percentage coverage of mpileup result/bowtie2inspect result *100
#Figures made in graphpad prism
#######################################################################################
#Calculate intersection for 3 K27ac BED files to make a venn diagram and find sites that are present in all AD samples
bedtools multiinter -header -i H3K27ac_4FC_AD3_DBsites_fixedmodel_11072018.bed H3K27ac_4FC_AD6_DBsites_fixedmodel_11072018.bed H3K27ac_4FC_AD12_DBsites_fixedmodel_11072018.bed -names K27ac_AD3 K27ac_AD6 K27ac_AD12 > K27ac_intersect.txt
#give a file that can be imported into R and then use gplots venn diagram
