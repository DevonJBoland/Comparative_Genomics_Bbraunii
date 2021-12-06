#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=MasucraMegaReads-Lrace        # job name
#SBATCH --time=02-00:00:00
#SBATCH --nodes=1             # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=80         # tasks (commands) per compute node
#SBATCH --cpus-per-task=1           # CPUs (threads) per command
#SBATCH --mem=2900G                    # total memory per node
#SBATCH --partition=bigmem          # use large (3TB) memory node
#SBATCH --output=MasucraMegaReads-Lrace.stdout.%j.txt          # save stdout to file
#SBATCH --error=MasucraMegaReads-Lrace.stderr.%j.txt           # save stderr to file

#SBATCH --mail-type=END
#SBATCH --mail-user=devonjboland@tamu.edu

##############
# Created on: 03-03-21
# Created by: Devon J. Boland
##############

# Load Modules
module load GCC/7.3.0-2.30 OpenMPI/3.1.1 MaSuRCA/4.0.1-Perl-5.28.0

# Run command from stdout of MASURCA-RUN.sh
sbatch -D /scratch/user/devonjboland/Canu-3TBNode-CanuAssembly/Masurca/mr_pass1 -J create_mega_reads -a 1-101 -n 80 -p bigmem -N 1 mr_pass1/create_mega_reads.sh
