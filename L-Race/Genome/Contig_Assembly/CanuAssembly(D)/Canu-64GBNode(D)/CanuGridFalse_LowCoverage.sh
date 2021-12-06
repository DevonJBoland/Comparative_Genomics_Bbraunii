#BSUB -L /bin/bash
#BSUB -J Canu_GridFalse_LowCoverage-Test3
#BSUB -o Canu_GridFalse_LowCoverage-Test3.txt
#BSUB -W 720:00
#BSUB -n 20
#BSUB -M 2700
#BSUB -R 'rusage[mem=2700]'
#BSUB -R 'span[ptile=20]'
#BSUB -e Canu_GridFalse_LowCoverage-Test3.stderr.txt

#BSUB -u devonjboland@tamu.edu
#BSUB -B -N


###############
# Created on: 01-29-2021
# Created by: Devon J. Boland
# Pass useGrid=remote in Canu executive wrapper script to generate a dry-run output
# to create job files to be sumbitted later with tamulauncher
###############

# Load Modules
module load Canu/2.0-intel-2017A-Perl-5.24.0

# Assign Variables and Directories
workdir=/scratch/user/devonjboland/CSH_ONT_ABLRaces_Bbraunii
cell1reads=$workdir/cell1/barcode10.fastq
cell2reads=$workdir/cell2/barcode10.fastq
mkdir $workdir/LowCoverage_CanuAssembly
outputdir=$workdir/LowCoverage_CanuAssembly

# Canu Command to Create Dryrun Scripts
canu useGrid=false maxMemory=54 maxThreads=20 minReadLength=600 rawErrorRate=0.500 correctedErrorRate=0.216 corMaxEvidenceErate=0.15 -d $outputdir -p LraceTest2 genomesize=211.3m -nanopore-raw $cell1reads $cell2reads
