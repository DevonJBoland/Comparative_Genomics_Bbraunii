#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=RNASeq-Trim-Map    #Set the job name to "JobExample2"
#SBATCH --time=02-00:00:00               #Set the wall clock limit to 6hr and 30min
#SBATCH --ntasks=80
#SBATCH --ntasks-per-node=80          #Request 20 tasks/cores per node
#SBATCH --mem=2900G                     #Request 2GB per node
#SBATCH --partition=bigmem
#SBATCH --output=stdout.%j.txt          # save stdout to file
#SBATCH --error=stderr.%j.txt           # save stderr to file 

##SET OPTIONAL DIRECTIVES
#SBATCH --mail-type=ALL              #Send email on all job events
#SBATCH --mail-user=devonjboland@tamu.edu    #Send all emails to


## Set Variables ##
#workdir='/scratch/user/devonjboland/L-Race/Transcriptome-Assembly'
#rnaseq='/scratch/user/devonjboland/L-Race/Transcriptome-Assembly/RNA-Seq'
genome='Coccomyxa_C169_v2_masked_genomic_scaffolds.fasta'
pre='FastQC-pre'
post='FastQC-post'
mkdir ${pre}
mkdir ${post}
memory=2900
threads=80

## Load Seqtk, Trimmomatic, FastQC Modules ##
ml GCC/7.3.0-2.30  OpenMPI/3.1.1 seqtk/1.3
module load Trimmomatic/0.39-Java-11
ml FastQC/0.11.9-Java-11
## Deinterleave RNA seq reads ##
#cd ${rnaseq}
#ls *.fastq > list.txt
cat list.txt | while read line
do
## Perform Initial FastQC Analysis
fastqc -o ${pre} -t ${threads} ${line}.fastq.gz
#seqtk seq $line -1 | pigz -p ${threads} -c > ${line}_R1.fq.gz
#seqtk seq $line -2 | pigz -p ${threads} -c > ${line}_R2.fq.gz
## Trim RNA-Seq Data Using Trimmomatic ##
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar SE -threads $threads ${line}.fastq.gz ${line}_trim ILLUMINACLIP:TruSeq.fa:2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:15
## FastQC Analysis of Trimmed RNA seq PE Reads ##
fastqc -o ${post} -t ${threads} ${line}_trim
done
ml purge

## Load Hisat2 Module(s) ##
ml GCC/7.3.0-2.30  OpenMPI/3.1.1 SAMtools/1.9 HISAT2/2.2.0
## Index Bbraunii Genome ##
#cd ${workdir}
hisat2-build --threads ${threads} ${genome} ht2
## Map Trimmed RNA Seq PE Reads to Bbraunii Genome ##
#cat ${rnaseq}/List.txt | while read line
#do
hisat2 -p ${threads} -x ht2 -U CG1_C_subellipsoidea_C-169_RNA-Seq_trim,CG2_C_subellipsoidea_C-169_RNA-Seq_trim,CG3_C_subellipsoidea_C-169_RNA-Seq_trim | samtools sort -o mapped-trim-reads.bam -T reads.tmp -
#done
#samtools merge mapped-trim-reads.bam B.braunii_race.L_Day.00.fastq_reads.sorted.bam B.braunii_race.L_Day.03.fastq_reads.sorted.bam B.braunii_race.L_Day.07.fastq_reads.sorted.bam B.braunii_race.L_Day.14.fastq_reads.sorted.bam B.braunii_race.L_Day.21.fastq_reads.sorted.bam B.braunii_race.L_Day.28.fastq_reads.sorted.bam
samtools index -@ $threads mapped-trim-reads.bam
ml purge
