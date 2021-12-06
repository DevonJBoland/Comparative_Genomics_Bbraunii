#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=BbrauniiLrace-Canu-3TBGrace        # job name
#SBATCH --time=2:00:00             # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=1           # CPUs (threads) per command
#SBATCH --mem=2G                    # total memory per node
#SBATCH --partition=bigmem          # use large (3TB) memory node
#SBATCH --output=BbrauniiLrace-Canu-3TBGrace.stdout.%j.txt          # save stdout to file
#SBATCH --error=BbrauniiLrace-Canu-3TBGrace.stderr.%j.txt           # save stderr to file

##############
# Created on: 02-05-21
# Created by: Devon J. Boland
##############

# Load Prerequesite Modules
module load GCCcore/8.3.0

# Load Canu Module
module load canu/2.1.1-Java-1.8

# Set Variables and Directories
workdir=/scratch/user/devonjboland/Canu-3TBNode-CanuAssembly
cell1reads=$workdir/ONTDataset/RaceL_cell1_barcode10.fastq
cell2reads=$workdir/ONTDataset/RaceL_cell2_barcode10.fastq
mkdir $workdir/Canu-3TBGraceNode-Output
outputdir=$workdir/Canu-3TBGraceNode-Output

#Execute Canu to Assemble Unitigs
canu useGrid=true gridOptions='--time=2-00:00:00 --partition=bigmem' rawErrorRate=0.500 correctedErrorRate=0.216 corMaxEvidenceErate=0.15 -d $outputdir -p LraceTest2 genomesize=211.3m -nanopore $cell1reads $cell2reads corMhapFilterThreshold=0.0000000002 corMhapOptions="--threshold 0.78 --num-hashes 768 --num-min-matches 2 --ordered-sketch-size 1000 --ordered-kmer-size 14 --min-olap-length 2000 --repeat-idf-scale 50"

# This job keeps exceeding my allocated disk usage space 5TB
# Searching the Canu Git page, I fdound a closed issue: https://github.com/marbl/canu/issues/587

# The problem seems to be caused by the high repeate nature of my genome
# To overcome I will try the the following code fix:
# corMhapFilterThreshold=0.0000000002 corMhapOptions="--threshold 0.78 --num-hashes 768 --num-min-matches 2 --ordered-sketch-size 1000 --ordered-kmer-size 14 --min-olap-length 2000 --repeat-idf-scale 50"