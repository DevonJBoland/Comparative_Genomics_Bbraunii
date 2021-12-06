#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=Trimmomatic       #Set the job name to "JobExample2"
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

# Load Modules
ml FastQC/0.11.9-Java-11
ml Trimmomatic/0.39-Java-11

# Trim Illumina Reads
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -t 80 Brace_CS358_R1.fastq.gz Brace_CS358_R2.fastq.gz -baseout Brace_CS358 ILLUMINACLIP:sequence-adapters.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
