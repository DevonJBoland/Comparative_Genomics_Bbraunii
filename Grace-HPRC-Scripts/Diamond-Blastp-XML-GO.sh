#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=Diamond-GO-Annotation    #Set the job name to "JobExample2"
#SBATCH --time=02-00:00:00           # max job run time dd-hh:mm:ss
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
#Arace='/scratch/user/devonjboland/GO-Annotation/BobraA.fasta'
#Lrace='/scratch/user/devonjboland/GO-Annotation/BobraL.fasta'
Brace='/scratch/user/devonjboland/GO-Annotation/Bbraunii_502_v2.1.protein.fa'
nr='/scratch/data/bio/diamond/nr.dmnd'
threads=$SLURM_CPUS_PER_TASK


# Load Blast Moduels
ml GCC/9.3.0 DIAMOND/2.0.11

# Use Diamond, instead of blast, it is much faster
#diamond blastp -q ${Arace} -d ${nr} -p ${threads} --fast -k 20 -f 5 -o Arace-Diamond.xml
#diamond blastp -q ${Lrace} -d ${nr} -p ${threads} --fast -k 20 -f 5 -o Lrace-Diamond.xml
diamond blastp -q ${Brace} -d ${nr} -p ${threads} --fast -k 20 -f 5 -o Brace-Diamond.xml