#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=BRAKER-UTR    #Set the job name to "JobExample2"
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

# Load Beast Module
ml GCC/7.3.0-2.30  OpenMPI/3.1.1  Beast/2.5.1

# Assign Variables
xml='Chlorophyta.xml'
threads=$SLURM_CPUS_PER_TASK

# Run Beast
java -Xms2000g -jar /sw/eb/sw/Beast/2.5.1-foss-2018b/lib/beast.jar -threads ${threads} -beagle_CPU -beagle_double ${xml}
