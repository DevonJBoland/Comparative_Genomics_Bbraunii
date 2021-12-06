#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=LRace-Illumina_AbyssAssembly        # job name
#SBATCH --time=2-00:00:00           # max job run time dd-hh:mm:ss
#SBATCH --nodes=1                   # total compute nodes
#SBATCH --ntasks-per-node=80        # tasks (commands) per compute node
#SBATCH --cpus-per-task=1           # CPUs (threads) per command
#SBATCH --mem=2900G                  # total memory per node
#SBATCH --partition=bigmem          # use large (3TB) memory node
#SBATCH --output=LRace-Illumina_AbyssAssembly.stdout.%j.txt          # save stdout to file
#SBATCH --error=LRace-Illumina_AbyssAssembly.stderr.%j.txt           # save stderr to file

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

################################### VARIABLES ##################################
# TODO Edit these variables as needed:
########## INPUTS ##########
pe1_1='/scratch/user/devonjboland/Canu-3TBNode-CanuAssembly/IlluminaDataSet/Race_L_R1.fastq.gz'
pe1_2='/scratch/user/devonjboland/Canu-3TBNode-CanuAssembly/IlluminaDataSet/Race_L_R2.fastq.gz'

######## PARAMETERS ########
mpi_threads=80              # value of (--nodes * --ntasks-per-node)
serial_threads=80            # value of --cpus-per-task
#kmer=111                    # max kmer available is 128
#min_pairs4contigs=5     # indicates 5 mate pairs needed to join contigs

########## OUTPUTS #########
prefix='Lrace_ILLUMINA'

################################### COMMANDS ###################################
# NOTE: running on too many nodes can actually result in a slower assembly
#gzip Race_L_R1.fastq
#gzip Race_L_R2.fastq

for k in `seq 50 8 90`; do
mkdir k$k
abyss-pe j=$serial_threads np=$mpi_threads -C k$k k=$k name=$prefix lib='lib1' lib1="$pe1_1 $pe1_2"
done
#abyss-fac k*/Lrace-contigs.fa
################################################################################
<<CITATION
    - Acknowledge TAMU HPRC: https://hprc.tamu.edu/research/citations.html

    - ABySS:
        Simpson, J. T., Wong, K., Jackman, S. D., Schein, J. E., Jones, S. J., & Birol, I. (2009).
        ABySS: a parallel assembler for short read sequence data. Genome research, 19(6), 1117-1123.
CITATION
