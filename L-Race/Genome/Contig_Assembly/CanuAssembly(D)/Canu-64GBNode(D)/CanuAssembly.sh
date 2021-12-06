#BSUB -L /bin/bash
#BSUB -J Canu_AssemblyCell2-Test2
#BSUB -o Canu_AssemblyCell2-Test2.txt
#BSUB -W 48:00
#BSUB -n 20
#BSUB -M 2700
#BSUB -R 'rusage[mem=2700]'
#BSUB -R 'span[ptile=20]'
#BSUB -e Canu_AssemblyCell2-Test2.stderr.txt

#BSUB -u devonjboland@tamu.edu
#BSUB -B -N


###############
# Created on: 01-15-2021
# Created by: Devon J. Boland
# Will only assemble cell1 and cell2 barcode10 L race reads
###############

# Load Modules
module load Canu/2.0-intel-2017A-Perl-5.24.0

# Assign Directories and Variables
workdir=/scratch/user/devonjboland/CSH_ONT_ABLRaces_Bbraunii
mkdir $workdir/Canu-Cell1Cell2-Barcode10-Lrace
outputdir=$workdir/Canu-Cell1Cell2-Barcode10-Lrace

# Assemble Fastq Reads using Canu
canu useGrid=false rawErrorRate=0.500 correctedErrorRate=0.144 -d $outputdir -p LraceTest1 genomesize=211.3m -nanopore $workdir/cell1/barcode10.fastq $workdir/cell2/barcode10.fastq