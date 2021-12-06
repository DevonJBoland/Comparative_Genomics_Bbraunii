#BSUB -L /bin/bash
#BSUB -J Canu_GridRemote-Test2
#BSUB -o Canu_GridRemote-Test2.txt
#BSUB -W 720:00
#BSUB -n 20
#BSUB -M 2700
#BSUB -R 'rusage[mem=2700]'
#BSUB -R 'span[ptile=20]'
#BSUB -e Canu_GridRemote-Test2.stderr.txt

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
outputdir=$workdir/DryRunScripts

# Canu Command to Create Dryrun Scripts
canu useGrid=false maxMemory=54 maxThreads=20 rawErrorRate=0.500 correctedErrorRate=0.144 -d $outputdir -p LraceTest2 genomesize=211.3m -nanopore $cell1reads $cell2reads
