#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=Mafft-RepeatLibrary-Arace   #Set the job name to "JobExample2"
#SBATCH --time=02-00:00:00               #Set the wall clock limit to 6hr and 30min
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=80          #Request 20 tasks/cores per node
#SBATCH --mem=2900G                     #Request 2GB per node
#SBATCH --partition=bigmem
#SBATCH --output=stdout.%j.txt          # save stdout to file
#SBATCH --error=stderr.%j.txt           # save stderr to file 

##SET OPTIONAL DIRECTIVES
#SBATCH --mail-type=ALL              #Send email on all job events
#SBATCH --mail-user=devonjboland@tamu.edu    #Send all emails to

# Variables to set
fasta='Filtered-consensi.fa'
msa='consensi-algn.fasta'
 

# Load ClustalOmega Moules & Align Initial Repeats to generate a phylip file
#ml iccifort/2020.1.217  impi/2019.7.217  MAFFT/7.471-with-extensions
#mafft --thread $SLURM_CPUS_PER_TASK ${fasta} > ${msa}
#ml purge

# Load RAxML to perform ML calculations
ml  iccifort/2019.5.281  impi/2018.5.288  RAxML/8.2.12-hybrid-avx2
raxmlHPC -f a -N 20 -m GTRGAMMA -p 12345 -x 12345 -s ${msa} -n Arace-InitialRepeat.tree
