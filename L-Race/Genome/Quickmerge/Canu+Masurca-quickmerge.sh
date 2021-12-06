#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=QuickmergeTest3       # job name
#SBATCH --time=02-00:00:00
#SBATCH --nodes=1             # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=80         # tasks (commands) per compute node
#SBATCH --cpus-per-task=1           # CPUs (threads) per command
#SBATCH --mem=2900G                    # total memory per node
#SBATCH --partition=bigmem          # use large (3TB) memory node
#SBATCH --output=stdout.%j.txt          # save stdout to file
#SBATCH --error=stderr.%j.txt           # save stderr to file

#SBATCH --mail-type=END
#SBATCH --mail-user=devonjboland@tamu.edu

##############
# Created on: 03-22-21
# Created by: Devon J. Boland
##############

# Load Modueles 
ml GCC/9.3.0 quickmerge/0.3

# Path to config files
Canu=/scratch/user/devonjboland/Canu-3TBNode-CanuAssembly/Assemblies/CanuLrace-PILON.contigs.fasta
Masurca=/scratch/user/devonjboland/Canu-3TBNode-CanuAssembly/Assemblies/MasurcaLRace-PILON.scaffolds.fasta
prefix=Merged-Canu-Lrace-scaffolds.

# Run quickmerge wrapper script
merge_wrapper.py -pre $prefix.ml5000. -l 1202580 -ml 5000 $Canu $Masurca
wait
merge_wrapper.py -pre $prefix.ml50000. -l 1202580 -ml 50000 $Canu $Masurca
wait
merge_wrapper.py -pre $prefix.ml500000. -l 1202580 -ml 500000 $Canu $Masurca
