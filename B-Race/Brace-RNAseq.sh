#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=Brace-RNAseq    #Set the job name to "JobExample2"
#SBATCH --time=00-10:00:00           # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=80          # CPUs (threads) per command
#SBATCH --mem=2900G                  # total memory per node
#SBATCH --partition=bigmem
#SBATCH --output=stdout.%j.txt          # save stdout to file
#SBATCH --error=stderr.%j.txt           # save stderr to file

##SET OPTIONAL DIRECTIVES
#SBATCH --mail-type=ALL              #Send email on all job events
#SBATCH --mail-user=devonjboland@tamu.edu    #Send all emails to

# Assign Variables
workdir='/scratch/user/devonjboland/B-Race'
rnaseq=${workdir}/'rnaseq'
truseq=${workdir}/'TruSeq.fa'
memory=2900
threads=$SLURM_CPUS_PER_TASK
pre=${workdir}/'pre-trim'
post=${workdir}/'post-trim'
genome=${workdir}/'Bbraunii_502_2.0.fa'

#########################################################################################
# Quality Analysis/Assesment Before Read Trimming # 
#########################################################################################
# Load FastQC modules
ml FastQC/0.11.9-Java-11
mkdir ${pre}
cat data.txt | while read line
do
fastqc -o ${pre} -t ${threads} ${rnaseq}/${line}
done
ml purge
#########################################################################################

#########################################################################################
# Trim Low Quality Bases and Adapters from Filtered Read Libraries # 
#########################################################################################
# Load Trimmomatic 
ml Trimmomatic/0.39-Java-11
# Trim off any adapter seqeunces from the Illumina library kits
cat data.txt | while read line
do
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar SE -threads $threads ${rnaseq}/${line} ${rnaseq}/${line}-trim ILLUMINACLIP:${truseq}:2:30:10
done
ml purge
#########################################################################################

#########################################################################################
# Quality Analysis/Assesment After Read Trimming # 
#########################################################################################
# Load FastQC modules
ml FastQC/0.11.9-Java-11
mkdir ${post}
ls ${rnaseq}/*-trim.fastq > trim.txt
cat trim.txt | while read line
do
fastqc -o ${post} -t ${threads} ${rnaseq}/${line}
done
ml purge
#########################################################################################

#########################################################################################
# Map Data Sets to B race Genome # 
#########################################################################################
## Load Hisat22 Module(s) ##
ml GCC/7.3.0-2.30  OpenMPI/3.1.1 SAMtools/1.9 HISAT2/2.2.0
## Index Bbraunii Genome ##
#hisat2-build --threads ${threads} ${genome} genome
## Map Trimmed RNA Seq PE Reads to Bbraunii Genome ##
cat trim.txt | while read line
do
hisat2 -p ${threads} -x genome -U ${line} | samtools sort -o ${line}_reads.sorted-masked.bam -T reads.tmp -
#done
cd ${rnaseq}
samtools merge Brace-mapped-trim-reads-masked.bam B.braunii_race.B_Day.00-trim.fastq_reads.sorted-masked.bam B.braunii_race.B_Day.03-trim.fastq_reads.sorted-masked.bam B.braunii_race.B_Day.07-trim.fastq_reads.sorted-masked.bam B.braunii_race.B_Day.14-trim.fastq_reads.sorted-masked.bam B.braunii_race.B_Day.21-trim.fastq_reads.sorted-masked.bam B.braunii_race.B_Day.28-trim.fastq_reads.sorted-masked.bam
samtools index -@ $threads Brace-mapped-trim-reads-masked.bam
ml purge
#########################################################################################