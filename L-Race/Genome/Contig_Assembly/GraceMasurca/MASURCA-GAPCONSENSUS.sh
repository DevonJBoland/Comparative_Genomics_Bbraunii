#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=MasucraGapConsensus-Lrace        # job name
#SBATCH --time=02-00:00:00
#SBATCH --nodes=1             # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=80         # tasks (commands) per compute node
#SBATCH --cpus-per-task=1           # CPUs (threads) per command
#SBATCH --mem=2900G                    # total memory per node
#SBATCH --partition=bigmem          # use large (3TB) memory node
#SBATCH --output=MasucraGapConsensus-Lrace.stdout.%j.txt          # save stdout to file
#SBATCH --error=MasucraConsensus-Lrace.stderr.%j.txt           # save stderr to file

#SBATCH --mail-type=END
#SBATCH --mail-user=devonjboland@tamu.edu

##############
# Created on: 03-03-21
# Created by: Devon J. Boland
##############

# Load Modules
module load GCC/7.3.0-2.30 OpenMPI/3.1.1 MaSuRCA/4.0.1-Perl-5.28.0

# Run Command in MASURCA stdout from assemble.sh
(cd mr.41.17.12.0.02.join_consensus.tmp && sbatch -D /scratch/user/devonjboland/Canu-3TBNode-CanuAssembly/Masurca/mr.41.17.12.0.02.join_consensus.tmp -J join_mega_reads -a 1-7 -n 80 -p bigmem -N 1 do_consensus.sh);
