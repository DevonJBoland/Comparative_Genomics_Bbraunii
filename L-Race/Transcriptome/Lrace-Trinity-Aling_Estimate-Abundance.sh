#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=Trinity-Aling-Estimate-Abundance    #Set the job name to "JobExample2"
#SBATCH --time=2-0:00:00               #Set the wall clock limit to 6hr and 30min
#SBATCH --ntasks=80
#SBATCH --ntasks-per-node=80          #Request 20 tasks/cores per node
#SBATCH --mem=2900G                     #Request 2GB per node
#SBATCH --partition=bigmem
#SBATCH --output=stdout.%j.txt          # save stdout to file
#SBATCH --error=stderr.%j.txt           # save stderr to file 

##SET OPTIONAL DIRECTIVES
#SBATCH --mail-type=ALL              #Send email on all job events
#SBATCH --mail-user=devonjboland@tamu.edu    #Send all emails to

## Variables ##
transcriptome='/scratch/user/devonjboland/L-Race/Transcriptome-Assembly/Lrace-Transcriptome_nofilter.fasta'
pe1-trim='/scratch/user/devonjboland/L-Race/Transcriptome-Assembly/Trinity-Assembly/RNA-Seq/L_race_RNA_1P.fastq.gz'
pe2-trim='/scratch/user/devonjboland/L-Race/Transcriptome-Assembly/Trinity-Assembly/RNA-Seq/L_race_RNA_2P.fastq.gz'
cat /scratch/user/devonjboland/L-Race/Transcriptome-Assembly/Trinity-Assembly/RNA-Seq/*U > /scratch/user/devonjboland/L-Race/Transcriptome-Assembly/Trinity-Assembly/RNA-Seq/L_race_RNA_U.fq
trimmed-unpaired='/scratch/user/devonjboland/L-Race/Transcriptome-Assembly/Trinity-Assembly/RNA-Seq/L_race_RNA_U.fq'
pe1='/scratch/user/devonjboland/L-Race/Transcriptome-Assembly/Trinity-Assembly/RNA-Seq/L_race_RNA_R1.fq.gz'
pe2='/scratch/user/devonjboland/L-Race/Transcriptome-Assembly/Trinity-Assembly/RNA-Seq/L_race_RNA_R2.fq.gz'
trim-paired='/scratch/user/devonjboland/L-Race/Transcriptome-Assembly/Abundance-Trinity/Trimmed-Paired/'
trim-unpaired='/scratch/user/devonjboland/L-Race/Transcriptome-Assembly/Abundance-Trinity/Trimmed-Unpaired/'
untrimmed='/scratch/user/devonjboland/L-Race/Transcriptome-Assembly/Abundance-Trinity/Untrimmed/'
memory=2900
threads=80

## Use Aling_and_Estimiate_abundance.pl/Trinity to Remove Low Expression Isoforms ##
## Load Trinity ##
ml GCC/8.3.0 OpenMPI/3.1.4 Trinity/2.10.0-Python-3.7.4
abundance='/sw/eb/sw/Trinity/2.10.0-foss-2019b-Python-3.7.4/trinityrnaseq-v2.10.0/util/align_and_estimate_abundance.pl'
## Run Estimation on Paired reads after trimming ##
perl ${abundance} --transcripts ${transcriptome} --seqType fq --left ${pe1} --right ${pe2} --est_method RSEM --aln_method bowtie2 --trinity_mode --prep_reference --output_dir ${trim-paired} --thread_count ${threads}
## Run Estimation on Unpaired (Single) reads after trimming ##
perl ${abundance} --transcripts ${transcriptome} --seqType fq --single ${trimmed-unpaired} --est_method RSEM --aln_method bowtie2 --trinity_mode --prep_reference --output_dir ${trim-unpaired} --thread_count ${threads}
## Run Estimation on Untrimmed reads ##
perl ${abundance} --transcripts ${transcriptome} --seqType fq --left ${pe1} --right ${pe2} --est_method RSEM --aln_method bowtie2 --trinity_mode --prep_reference --output_dir ${untrimmed} --thread_count ${threads}
