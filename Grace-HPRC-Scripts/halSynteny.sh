#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=hal2synteny-Test    #Set the job name to "JobExample2"
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

# Load HAL toolkit
ml GCC/10.2.0  OpenMPI/4.0.5 HAL/2.1

# Assign Variables
#workdir='/scratch/user/devonjboland/Cactus/halSynteny'
Arace='BobraA'
Brace='Botrbrau1'
Lrace='BobraL'
hal='/scratch/user/devonjboland/Cactus/Braunii-Cactus/algae.hal'

# Create PSL files for B race to L race alinments
hal2synteny --queryGenome ${Lrace} --targetGenome ${Brace} --maxAnchorDistance 5000 --minBlockSize 5000 ${hal} BvsL2.psl
# Create PSL files for B race to A race alignments
hal2synteny --queryGenome ${Arace} --targetGenome ${Brace} --maxAnchorDistance 5000 --minBlockSize 5000 ${hal} BvsA2.psl
# Create PSL files for L race to A race alignments
hal2synteny --queryGenome ${Arace} --targetGenome ${Lrace} --maxAnchorDistance 5000 --minBlockSize 5000 ${hal} LvsA2.psl
