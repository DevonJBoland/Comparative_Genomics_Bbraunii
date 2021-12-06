#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=Cactus-Bbraunii            # job name
#SBATCH --time=2-00:00:00           # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=80          # CPUs (threads) per command
#SBATCH --mem=2900G                  # total memory per node
#SBATCH --partition=bigmem
#SBATCH --output=stdout.%j.txt          # save stdout to file
#SBATCH --error=stderr.%j.txt           # save stderr to file

##SET OPTIONAL DIRECTIVES
#SBATCH --mail-type=ALL              #Send email on all job events
#SBATCH --mail-user=devonjboland@tapestatmu.edu    #Send all emails to

jobStorePath='/scratch/user/devonjboland/Cactus/jobStore'
seqFile='/scratch/user/devonjboland/Cactus/seqFile'
output='/scratch/user/devonjboland/Cactus/Cactus.hal'
threads=$SLURM_CPUS_PER_TASK
memory=2900G

# Execute Cactus via Singularity Container
export SINGULARITY_BINDPATH="/scratch,$TMPDIR"
singularity exec /sw/hprc/sw/bio/containers/cactus-v2.0.3.sif cactus ${jobStorePath} ${seqFile} ${output}
