#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=canu_grid        # job name
#SBATCH --time=01:00:00             # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=1           # CPUs (threads) per command
#SBATCH --mem=2G                    # total memory per node
#SBATCH --partition=bigmem          # use large (3TB) memory node
#SBATCH --output=stdout.%j          # save stdout to file
#SBATCH --error=stderr.%j           # save stderr to file

module load GCCcore/8.3.0
module load canu/2.1.1-Java-1.8

<<README
    - Canu Tutorial: http://canu.readthedocs.io/en/latest/tutorial.html
README

################################### VARIABLES ##################################
# TODO Edit these variables as needed:
########## INPUTS ##########
pacbio_raw_reads='/scratch/data/bio/GCATemplates/pacbio/neurospora_crassa/OR74A_filtered_subreads.fastq'

######## PARAMETERS ########
genome_size='40m'             	    # supported units: g, m, k
stop_on_low_coverage='stopOnLowCoverage=4'  # default 10, using 4 for sample dataset, adjust as needed

########## OUTPUTS #########
prefix='n_crassa'
assembly_directory='grid_build_bigmem_2.1.1_out'

################################### COMMANDS ###################################
# bigmem partition time limit is 2 days
canu useGrid=true gridOptions='--time=2-00:00:00 --partition=bigmem' \
-p $prefix -d $assembly_directory genomeSize=$genome_size \
 -pacbio-raw $pacbio_raw_reads $stop_on_low_coverage



<<CITATION
    - Acknowledge TAMU HPRC: https://hprc.tamu.edu/research/citations.html

    - Canu: Koren S, Walenz BP, Berlin K, Miller JR, Phillippy AM.
            Canu: scalable and accurate long-read assembly via adaptive
            k-mer weighting and repeat separation. Genome Research. (2017).
CITATION
