#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=OrthoFinder-Algae    #Set the job name to "JobExample2"
#SBATCH --time=02-00:00:00               #Set the wall clock limit to 6hr and 30min
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=80          #Request 20 tasks/cores per node
#SBATCH --mem=2900G                     #Request 2GB per node
#SBATCH --partition=bigmem
#SBATCH --output=stdout.%j.txt          # save stdout to file
#SBATCH --error=stderr.%j.txt           # save stderr to file 

##SET OPTIONAL DIRECTIVES
#SBATCH --mail-type=ALL              #Send email on all job events
#SBATCH --mail-user=devonjboland@tamu.edu    #Send all emails to

# Assign Variables
proteins='/scratch/user/devonjboland/OrthoFinder/Proteins'
threads=$SLURM_CPUS_PER_TASK
# Load OrthoFinder
ml iccifort/2019.5.281  impi/2018.5.288 OrthoFinder/2.3.11-Python-3.7.4
# Run OrthoFinder
orthofinder -t ${threads} -a 20 -f ${proteins}
