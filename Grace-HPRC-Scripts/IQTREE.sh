#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=Otholog-BUSCO-Aln-Phylogenomic-Analysis            # job name
#SBATCH --time=2-00:00:00           # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=16         # tasks (commands) per compute node
#SBATCH --cpus-per-task=5          # CPUs (threads) per command
#SBATCH --mem=2900G                  # total memory per node
#SBATCH --partition=bigmem
#SBATCH --output=stdout.%j.txt          # save stdout to file
#SBATCH --error=stderr.%j.txt           # save stderr to file

# Assign Variables
sequences='/scratch/user/devonjboland/Phytozome/Missing-No-More-Than-3'
mafft='/scratch/user/devonjboland/Phytozome/MAFFT'
iqtree='/scratch/user/devonjboland/Phytozome/IQ-TREE'


#Load Mafft Modules
ml iccifort/2020.1.217  impi/2019.7.217  MAFFT/7.471-with-extensions
# Perform alignments on extracted BUSCO protein sequences from algae
cat aln.txt | while read line
do
mafft --thread ${threads} ${sequences}/${line} > ${mafft}/aln-${line}
done
ml purge

# Load IQ-TREE Modules
ml GCC/9.3.0  OpenMPI/4.0.3   IQ-TREE/2.1.2
# Perform ML-Tree Analysis
tamulauncher --commands-pernode 16 commands.txt
