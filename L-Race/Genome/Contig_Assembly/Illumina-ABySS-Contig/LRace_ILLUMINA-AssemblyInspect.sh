#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=LRace-Illumina_AbyssAssembly_Inspect        # job name
#SBATCH --time=2-00:00:00           # max job run time dd-hh:mm:ss
#SBATCH --nodes=1                   # total compute nodes
#SBATCH --ntasks-per-node=80        # tasks (commands) per compute node
#SBATCH --cpus-per-task=1           # CPUs (threads) per command
#SBATCH --mem=2900G                  # total memory per node
#SBATCH --partition=bigmem          # use large (3TB) memory node
#SBATCH --output=LRace-Illumina_AbyssAssemblyInspect.stdout.%j.txt          # save stdout to file
#SBATCH --error=LRace-Illumina_AbyssAssemblyInspect.stderr.%j.txt           # save stderr to file

#SBATCH --mail-type=END
#SBATCH --mail-user=devonjboland@tamu.edu

##############
# Created on: 02-25-21
# Created by: Devon J. Boland
##############

# Load Modules for ABySS
ml GCC/8.3.0
ml OpenMPI/3.1.4
ml ABySS/2.1.5

abyss-fac k*/Lrace_ILLUMINA-contigs.fa
