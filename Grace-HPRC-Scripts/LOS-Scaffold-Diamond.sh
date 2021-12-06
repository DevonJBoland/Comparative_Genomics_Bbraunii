#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=Diamond-LOS-Scaffold-Genes    #Set the job name to "JobExample2"
#SBATCH --time=00-01:00:00           # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=10          # CPUs (threads) per command
#SBATCH --mem=320G                  # total memory per node
#SBATCH --partition=medium
#SBATCH --output=stdout.%j.txt          # save stdout to file
#SBATCH --error=stderr.%j.txt           # save stderr to file

##SET OPTIONAL DIRECTIVES
#SBATCH --mail-type=ALL              #Send email on all job events
#SBATCH --mail-user=devonjboland@tamu.edu    #Send all emails to

# Assign Variables
query='LOS-Scaffold-Genes.faa'
nr='/scratch/data/bio/diamond/nr.dmnd'
threads=$SLURM_CPUS_PER_TASK

# Load Modules
ml GCC/9.3.0 DIAMOND/2.0.11

# Run Diamond Alignment against nt
diamond blastp -q ${query} -d ${nr} -p ${threads} --fast -k 5 -f 6 qseqid qlen sseqid slen qstart qend start send length evalue pident staxids sscinames stitle -o LOS-Scaffold-Genes-Diamond.txt
