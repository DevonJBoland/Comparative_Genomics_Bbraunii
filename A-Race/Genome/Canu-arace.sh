#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=Arace-Canu      #Set the job name to "JobExample2"
#SBATCH --time=2:00:00             # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=1           # CPUs (threads) per command
#SBATCH --mem=2G                    # total memory per node
#SBATCH --partition=bigmem          # use large (3TB) memory node
#SBATCH --output=stdout.%j.txt          # save stdout to file
#SBATCH --error=stderr.%j.txt           # save stderr to file

##SET OPTIONAL DIRECTIVES
#SBATCH --mail-type=ALL              #Send email on all job events
#SBATCH --mail-user=devonjboland@tamu.edu    #Send all emails to

## Set Variables and Directories ##
ontreads='/scratch/user/devonjboland/A-Race/trimmed-Arace-cell1-2-barcode08-ont.fastq.gz'
outputdir='/scratch/user/devonjboland/A-Race/Canu/'

## Load Modules ##
module load GCCcore/8.3.0
module load canu/2.1.1-Java-1.8
## Execute Canu to Assemble Unitigs ##
canu useGrid=true gridOptions='--time=2-00:00:00 --partition=bigmem' rawErrorRate=0.500 correctedErrorRate=0.216 corMaxEvidenceErate=0.15 -d $outputdir -p Brace-Canu genomesize=166m -nanopore $ontreads

