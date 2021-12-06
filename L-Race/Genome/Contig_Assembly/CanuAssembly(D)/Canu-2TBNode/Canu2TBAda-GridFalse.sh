#BSUB -L /bin/bash
#BSUB -J Canu-2TBNode-GridFalse
#BSUB -o Canu-2TBNode-GridFalse.%J.stdout.txt
#BSUB -e Canu-2TBNode-GridFalse.%J.stderr.txt
#BSUB -n 40
#BSUB -R "select[mem2tb]"
#BSUB -R "span[ptile=40]"
#BSUB -R "rusage[mem=49750]"
#BSUB -M 49750
#BSUB -q xlarge
#BSUB -W 240:00

#############################
# Created on: 02-04-2021
# Created by: Devon J. Boland
#############################

# Load Westmere module for Large >1TB memory jobs
module load Westmere
module load Canu/2.0-intel-2017A-Perl-5.24.0

# Set Varaibles and Make Directories
workdir=/scratch/user/devonjboland/CSH_ONT_ABLRaces_Bbraunii
cell1reads=$workdir/cell1/barcode10.fastq
cell2reads=$workdir/cell2/barcode10.fastq
mkdir $workdir/LowCoverage_CanuAssembly
outputdir=$workdir/LowCoverage_CanuAssembly

# Initiate Canu Executive Wrapper
canu useGrid=remote maxMemory=1990 maxThreads=40 rawErrorRate=0.500 correctedErrorRate=0.216 corMaxEvidenceErate=0.15 corMhapFilterThreshold=0.0000000002 corMhapOptions="--threshold 0.78 --num-hashes 768 --num-min-matches 2 --ordered-sketch-size 1000 --ordered-kmer-size 14 --min-olap-length 2000 --repeat-idf-scale 50" -d $outputdir -p LraceTest2 genomesize=211.3m -nanopore $cell1reads $cell2reads
