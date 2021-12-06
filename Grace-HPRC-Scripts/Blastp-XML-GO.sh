#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=blastp-GO-Annotation    #Set the job name to "JobExample2"
#SBATCH --time=02-00:00:00           # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=5         # tasks (commands) per compute node
#SBATCH --cpus-per-task=16          # CPUs (threads) per command
#SBATCH --mem=2900G                  # total memory per node
#SBATCH --partition=bigmem
#SBATCH --output=stdout.%j.txt          # save stdout to file
#SBATCH --error=stderr.%j.txt           # save stderr to file

##SET OPTIONAL DIRECTIVES
#SBATCH --mail-type=ALL              #Send email on all job events
#SBATCH --mail-user=devonjboland@tamu.edu    #Send all emails to

# Assign Variables
Arace='/scratch/user/devonjboland/GO-Annotation/Arace'
Lrace='/scratch/user/devonjboland/GO-Annotation/Lrace'
nr='/scratch/data/bio/blastdb-2021.05.08/nr'
threads=$SLURM_CPUS_PER_TASK


# Load Blast Moduels
ml GCC/10.2.0  OpenMPI/4.0.5  BLAST+/2.11.0

tamulauncher --commands-pernode 5 commands.txt
