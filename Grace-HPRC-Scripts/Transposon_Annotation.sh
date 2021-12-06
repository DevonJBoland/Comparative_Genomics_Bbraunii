#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=RepeatAnnotation            # job name
#SBATCH --time=2-00:00:00           # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=80          # CPUs (threads) per command
#SBATCH --mem=2900G                  # total memory per node
#SBATCH --partition=bigmem
#SBATCH --output=stdout.%j.txt          # save stdout to file
#SBATCH --error=stderr.%j.txt           # save stderr to file

# Set Variables
arace='Arace-Filtered-Transposon_Models_v1.fa'
table='hmmscan-tbl'
txt='hmmscan-txt'
lrace='Lrace-Filtered-Transposon-Models_v1.fa'
threads=$SLURM_CPUS_PER_TASK

# Load HMMER Modules
ml  GCC/10.2.0  OpenMPI/4.0.5  HMMER/3.3.2

# Perform HMMER Analysis for closest homolog in Dfam database
#A Race Filtered Repeat Library
hmmscan --cpu ${threads} --tblout ${table}-Arace --noali Dfam_curatedonly.hmm ${arace} 

# L Race Filtered Repeat Library
hmmscan --cpu ${threads} --tblout ${table}-Lrace --noali Dfam_curatedonly.hmm ${lrace}
