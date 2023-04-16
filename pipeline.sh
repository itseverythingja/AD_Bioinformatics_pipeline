#!/bin/bash
# In order for us to execute our NGS pipeline to perform read alignment, variant discovery and annotation, we started
# by organising our files. First, we made a directory for the DNA sequencing within the folder 'advanced_bio_assignment'
# using the mkdir command
mkdir advanced_bio_assignment
mkdir advanced_bio_assignment/dnaseq

# Next, we made the following studcture within the directory to keep the files organised
cd ~/advanced_bio_assignment/dnaseq
mkdir data raw results commands

# Within the folder 'data', we created 2 subdirectories: 'untrimmed_fastq' and 'trimmed_fastq'
cd ~/advanced_bio_assignment/dnaseq/data
mkdir untrimmed_fastq
mkdir trimmed_fastq

# The data, provided in the instructions, was the downloaded using the 'wget' command line
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R1.fastq.qz
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R2.fastq.qz
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/annotation.bed

# The raw fastq files were the coppied into the 'untrimmed_fastq' file and the bed file in 'data' directory'
mv *fastq.gz ~/advanced_bio_assignment/dnaseq/data/untrimmed_fastq
mv annotation.bed ~/advanced_bio_assignment/dnaseq/data

# Since we will be performing alignement later during the workflow, we downloaded the reference file
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz

# And saved the file in the 'data' directory
mv hg19.fa.gz ~/ngs_course/dnaseq/data/


# Before starting the pipeline, the tools required to for the analysis were downloaded
# First Anaconda was installed
wget https://repo.anaconda.com/archive/Anaconda3-2022.10-Linux-x86_64.sh

# In order for to execute Anaconda, we need to set the executable bit for the scridpt. This is a necessary step
chmod +x ./Anaconda3-2022.10-Linux-x86_64.sh

# In order to run the Anaconda installer script on a Linux operating system using the Bash shell, the following command
# line is used
bash ./Anaconda3-2022.10-Linux-x86_64.sh

# We reload the current user's Bash shell configuration file, we use this command line. This is useful since we  have
#  made changes to the configuration file and want to apply them immediately without having to log out and log back in 
source ~/.bashrc

# We add the 'defaults' 'biconda' and 'conda-forge' channels to the list of channels that conda searches for packages
# using the 'conda config' command line: used to configure various settings for conda
# '--add channels' option tells Conda to add a new channel to its list of channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# And lastly, we download all the tools needed for analysis on conda
conda install samtools
conda install bwa
conda install freebayes
conda install picard
conda install bedtools
conda install trimmomatic
conda install fastqc
conda install vcflib

# Now that the files are organised, the data required downloaded and stored, and the tools are 
# downloaded, we start with the quality assessment. This is done using the command line 'fastqc'
# which allows us to do come quality control check on raw sequence data
# first, we make sure we are in the right directory
cd ~/advanced_bio_assignment/dnaseq/data/untrimmed_fastq
# before running fastqc, we need to compress the files (.qz to .gz)
gzip NGS0001.R1.fastq.qz
gunzip -c NGS0001.R1.fastq.qz.gz > NGS0001.R1.fastq.gz
gzip NGS0001.R2.fastq.qz
gunzip -c NGS0001.R2.fastq.qz.gz > NGS0001.R2.fastq.gz
# then we run fastqc command line using multiple file names
fastqc *.fastq.gz

# We then create a home for our results 
mkdir ~/advanced_bio_assignment/dnaseq/results/fastqc_untrimmed_reads
# and move the results there 
mv *fastqc* ~/advanced_bio_assignment/dnaseq/results/fastqc_untrimmed_reads/

# To look at the results, we use FileZila, a seperate software downloaded on the computer to view
# the html files generated. However, the files need to be unpacked using 'unzip' program
# So first, we get to the right directory
cd ~/advanced_bio_assignment/dnaseq/results/fastqc_untrimmed_reads/
# and then we unpack using a 'for' loop to iterate through the list of the files in '*.zip'
$ for zip in *.zip
> do
> unzip $zip
> done

# to check the information contained in the unzipped files, we use these command lines
# for the first file
ls -lh NGS0001.R1_fastqc
# for the second file 
ls -lh NGS0001.R2_fastqc
# followed by 
head NGS0001.R1_fastqc/summary.txt
head NGS0001.R2_fastqc/summary.txt
# and we save all the texts into a 'full_report.txt' and move it to the 'raw' directory
cat */summary.txt > ~/advanced_bio_assignment/dnaseq/raw/fastqc_summaries.txt

## now we move onto trimming where we improve the quality of our reads, by trimming off any "bad" bases.
# we will use Trimmomatic to trim away adapters and filter out poor quality score reads
# trimmomatic is Java based and these are its usage instructions

# the 'PE' argument is a keyword that specifies we are working with paired-end reads
# we have to specify the '-threads' parameter because Trimmomatic uses all threads on a node by default
# For our paired-end fastq files, the command is
trimmomatic PE \
-threads 4 \
-phred33 \
/home/ubuntu/advanced_bio_assignment/dnaseq/data/untrimmed_fastq/NGS0001.R1.fastq.gz /home/ubuntu/advanced_bio_assignment/dnaseq/data/untrimmed_fastq/NGS0001.R2.fastq.gz \
-baseout /home/ubuntu/advanced_bio_assignment/dnaseq/data/trimmed_fastq/NGS0001.trimmed.R \
ILLUMINACLIP:/home/ubuntu/anaconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa:2:30:10 \
TRAILING:25 MINLEN:50

## now we assess the quality again
# but first we get back to the right directory
cd ~/advanced_bio_assignment/dnaseq/data/trimmed_fastq/
# and then we run fastqc
fastqc NGS0001.trimmed.R_1P  NGS0001.trimmed.R_1U  NGS0001.trimmed.R_2P  NGS0001.trimmed.R_2U
# now we place these files in a new directory
mkdir ~/advanced_bio_assignment/dnaseq/results/fastqc_trimmed_reads
mv *fastqc* ~/advanced_bio_assignment/dnaseq/results/fastqc_trimmed_reads

# To look at the results, we use FileZila, a seperate software downloaded on the computer to view
# the html files generated. However, the files need to be unpacked using 'unzip' program
# So first, we get to the right directory
cd ~/advanced_bio_assignment/dnaseq/results/fastqc_trimmed_reads
# and then we unpack using a 'for' loop to iterate through the list of the files in '*.zip'
$ for zip in *.zip
> do
> unzip $zip
> done

# to check the information contained in the unzipped files of the paired end, we use these command lines
ls -lh NGS0001.trimmed.R_1P_fastqc NGS0001.trimmed.R_2P_fastqc
# followed by
head NGS0001.trimmed.R_1P_fastqc/summary.txt
head NGS0001.trimmed.R_2P_fastqc/summary.txt
# and we save all the texts into a 'full_report.txt' and move it to the 'raw' directory
cat */summary.txt > ~/advanced_bio_assignment/dnaseq/raw/fastqc_trimmed_summaries.txt

## Read Alignment
# Now that we have our quality-trimmed reads, we can move on to read alignment. We perform read alignment using BWA MEM
# Firts, we make the index
mkdir -p ~/advanced_bio_assignment/dnaseq/data/reference
mv ~/advanced_bio_assignment/dnaseq/data/hg19.fa.gz ~/advanced_bio_assignment/dnaseq/data/reference/
bwa index ~/advanced_bio_assignment/dnaseq/data/reference//hg19.fa.gz
# Files generated by BWA index 
ls ~/advanced_bio_assignment/dnaseq/data/reference
hg19.fa.gz hg19.fa.gz.amb hg19.fa.gz.ann hg19.fa.gz.pac

## Read Alignment
# Now that we have our quality-trimmed reads, we can move on to read alignment. We perform read alignment using BWA MEM
# Firts, we make the index
mkdir -p ~/advanced_bio_assignment/dnaseq/data/reference
mv ~/advanced_bio_assignment/dnaseq/data/hg19.fa.gz ~/advanced_bio_assignment/dnaseq/data/reference/
bwa index ~/advanced_bio_assignment/dnaseq/data/reference//hg19.fa.gz
# Files generated by BWA index 
ls ~/advanced_bio_assignment/dnaseq/data/reference
hg19.fa.gz hg19.fa.gz.amb hg19.fa.gz.ann hg19.fa.gz.pac

# using the reference genome, run bwa mem 
bwa mem
-t 4
-v 1
-R '@RG\tID:HWI-D0011.50.H7AP8ADXX.1.WES01\tSM:WES01\tPL:ILLUMINA\tLB:nextera-wes01-blood\tDT:2017-02-23\tPU:HWI-D00119' -I 250,50
~/advanced_bio_assignment/dnaseq/data/reference/hg19.fa.gz
~/advanced_bio_assignment/dnaseq/data/trimmed_fastq/NGS0001.trimmed.R_1P
~/advanced_bio_assignment/dnaseq/data/trimmed_fastq/NGS0001.trimmed.R_2P >
~/advanced_bio_assignment/dnaseq/data/aligned_data/NGS0001.trimmed.R.sam

# Change directories into the aligned_data folder
cd ~/advanced_bio_assignment/dnaseq/data/aligned_data
# now we convert the sam file into bam format, sort it and generate an index using samtools
#IMPORTANT: depending on your samtools version, you might need to add the -S option to make samtools accept the sam input file 
$ samtools view -h -b NGS0001.trimmed.R.sam > NGS0001.trimmed.R.bam
   # with the above step, I kept getting this error message 'samtools: error while loading shared libraries: libcrypto.so.1.0.0: cannot open shared object file: No such file or directory'
   # so to fix it, the following steps are necessary
   # Run the following command to update your package sources
   sudo apt-get update
   # Finally, run the following command to install the package
   sudo apt-get install libssl1.0.0
   # that's when samtools command line (above) worked

#IMPORTANT: depending on your samtools version, you might need to add the -o option to make samtools direct the output to stdout
$ samtools sort NGS0001.trimmed.R.bam > NGS0001.trimmed.R.sorted.bam
# this will generate a .bai index file
$ samtools index NGS0001.trimmed.R.sorted.bam

$ ls
NGS0001.trimmed.R.bam  NGS0001.trimmed.R.sam  NGS0001.trimmed.R.sorted.bam  NGS0001.trimmed.R.sorted.bam.bai


## Post Alignment QC and Filtering
# We will Use 1) NGS: Picard > MarkDuplicates and then 2) SAMtools > Filter SAM or BAM, output SAM or BAM files to make a new BAM file which flags duplicate molecules and excludes poor quality reads
# MarkDuplicates
# We will use Picard to mark duplicated reads
picard MarkDuplicates I=NGS0001.trimmed.R.sorted.bam O=NGS0001.trimmed.R.sorted.marked.bam M=marked_dup_metrics.txt
samtools index NGS0001.trimmed.R.sorted.marked.bam

# Filter BAM based on mapping quality and bitwise flags using samtools
# We are going to filter the reads according to the following criteria
samtools view -F 1796  -q 20 -o NGS0001.trimmed.R.sorted.filtered.bam NGS0001.trimmed.R.sorted.marked.bam
samtools index NGS0001.trimmed.R.sorted.filtered.bam

## Alignment Statistics
# flagstat using samtools: generate flagstat statistics for the aligned BAM file
samtools flagstat NGS0001.trimmed.R.sorted.marked.bam
# to save the results in a file 
samtools flagstat NGS0001.trimmed.R.sorted.filtered.bam > NGS0001.flagstat.txt

# to generate a standard alignment statistics for a BAM file using
samtools idxstats NGS0001.trimmed.R.sorted.filtered.bam
# we save the file using the following command line
samtools idxstats NGS0001.trimmed.R.sorted.filtered.bam > NGS0001.idxstats.txt

# Determine the distribution of insert sizes using picard, 
picard CollectInsertSizeMetrics \
I=NGS0001.trimmed.R.sorted.filtered.bam \
O=NGS0001_size_metrics.txt \
H=NGS0001_size_histogram.pdf \
M=0.5
# determining the depth Coverage using bedtools, here the '-a' is an option that allows us to includ positions with zero coverage in the output, while omitting this option will only output
#  positions with non-zero coverage
samtools depth -a NGS0001.trimmed.R.sorted.filtered.bam > NGS0001.depth.txt


### Variant calling with Freebayes
#  We are going to convert to text format the reference 
zcat ~/advanced_bio_assignment/dnaseq/data/reference/hg19.fa.gz > ~/advanced_bio_assignment/dnaseq/data/reference/hg19.fa

# and now index it with samtools faidx
samtools faidx ~/advanced_bio_assignment/dnaseq/data/reference/hg19.fa

# calling variants with Freebayes
freebayes --bam ~/advanced_bio_assignment/dnaseq/data/aligned_data/NGS0001.trimmed.R.sorted.filtered.bam --fasta-reference ~/advanced_bio_assignment/dnaseq/data/reference/hg19.fa \
--vcf ~/advanced_bio_assignment/dnaseq/results/NGS0001.vcf

# compressing the resulting variant file (VCF)
bgzip ~/advanced_bio_assignment/dnaseq/results/NGS0001.vcf
   ## note that tabix had to be intall using the following command 
     sudo apt install tabix

# and finally index the VCF with tabix
tabix -p vcf ~/advanced_bio_assignment/dnaseq/results/NGS0001.vcf.gz


### Filtering the VCF
# We will apply the following suggested freebayes hard filter for human diploid sequencing and use vcffilter to apply it:
# QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1, where QUAL > 1: removes horrible sites QUAL / AO > 10 : additional contribution of each obs should be 10 log units
#  (~ Q10 per read) SAF > 0 & SAR > 0 : reads on both strands RPR > 1 & RPL > 1 : at least two reads “balanced” to each side of the site
vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" \
~/advanced_bio_assignment/dnaseq/results/NGS0001.vcf.gz > ~/advanced_bio_assignment/dnaseq/results/NGS0001_filtered.vcf

# using bedtools we can filter the vcf file for the regions
# the option intersect finds overlapping intervals in various ways, the option -wa writes the original entry in A for each overlap, 
bedtools intersect -header -wa -a ~/advanced_bio_assignment/dnaseq/results/NGS0001_filtered.vcf -b ~/advanced_bio_assignment/dnaseq/data/annotation.bed \
 ~/advanced_bio_assignment/dnaseq/results/NGS0001_filtered.vcf
bgzip ~/advanced_bio_assignment/dnaseq/results/NGS0001_filtered.vcf
tabix -p vcf ~/dvanced_bio_assignment/dnaseq/results/NGS0001_filtered.vcf.gz
