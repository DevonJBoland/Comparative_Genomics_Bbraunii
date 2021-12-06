#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=Quast-Analysis       # job name
#SBATCH --time=10:00:00             # max job run time dd-hh:mm:ss
#SBATCH --nodes=1
#SBATCH --ntasks=20       # tasks (commands) per compute node
#SBATCH --mem=360G                    # total memory per node
#SBATCH --output=stdout.%j.txt          # save stdout to file
#SBATCH --error=stderr.%j.txt           # save stderr to file

#SBATCH --mail-type=END
#SBATCH --mail-user=devonjboland@tamu.edu

##############
# Created on: 02-11-21
# Created by: Devon J. Boland
##############

# Load Modules
module load GCC/9.3.0
module load OpenMPI/4.0.3
module load QUAST/5.0.2-Python-3.8.2

# Set Variables
workdir=/scratch/user/devonjboland/Canu-3TBNode-CanuAssembly
outputdir=$workdir/Quast/Quickmerge
quast=/sw/eb/sw/QUAST/5.0.2-foss-2020a-Python-3.8.2/bin/quast

# Assemblies
Merged-5000=$workdir/Assemblies/Merged-Canu-Masurca-Lrace.ml5000.scaffolds.fasta
Merged-50000=$workdir/Assemblies/Merged-Canu-Masurca-Lrace.ml50000.scaffolds.fasta
Merged-500000=$workdir/Assemblies/Merged-Canu-Masurca-Lrace.ml500000.scaffolds.fasta
Merged-5000-Pilon=$workdir/Assemblies/
Merged-50000-Pilon=$workdir/Assemblies/
Merged-500000-Pilon=$workdir/Assemblies/

# Execute Quast
python $quast -t 20 -o $outputdir ${Merged-5000} ${Merged-50000} ${Merged-500000} ${Merged-5000-Pilon} ${Merged-50000-Pilon} ${Merged-500000-Pilon}
 
