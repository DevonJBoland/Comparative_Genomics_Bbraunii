#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=HAL-Analysis-Cactus            # job name
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

# Assign Variables
hal='algae.hal'
gffb='Botrbrau1_GeneCatalog_genes_20200805.gff'

# Load HAL/2.1
ml GCC/10.2.0  OpenMPI/4.0.5  HAL/2.1

# Validate HAL file
halValidate ${hal}
# Obtain HAL file statistics
halStats ${hal}
# Count Mutations
halSummarizeMutations ${hal}
ml purge

# Convert B race GFF file to BED using BEDOPS
# Load BEDOPS/2.4
ml GCC/7.3.0-2.30  OpenMPI/3.1.1  BEDOPS/2.4.35
# Convert
gff2bed < ${gffb} > ${gffb}.bed
ml purge

# Reload HAL/2.1
ml GCC/10.2.0  OpenMPI/4.0.5  HAL/2.1
# Calculate syntenny between B race Genes and A and L race genomes
# A Race
halLiftover ${hal} Botrbrau1 ${gffb}.bed BobraA BobraA_annotated.bed
# L Race
halLiftover ${hal} Botrbrau1 ${gffb}.bed BobraL BobraL_annotated.bed
