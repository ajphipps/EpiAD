#Download NGS data from AGRF through wget
#mount disk for datastorage
sudo mkdir /NGS_Data
sudo mount /dev/vdb /NGS_Data -t auto
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
ls 154467******K4me3*.fastq
cat 154467******K4me3*.fastq > 154467K4me3.fastq
#fastqc each resulting sample and record QC data
fastqc *all.fastq -o /NGS_Data/Andrew/CAGRF13394_QC
#also fastqc lane data for future reference and look for lane variation between samples
mkdir /NGS_Data/Andrew/CAGRF13394_QC/laneqc/
fastqc *R1.fastq -o /NGS_Data/Andrew/CAGRF13394_QC/laneqc/
