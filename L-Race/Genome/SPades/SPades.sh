#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=SPadesHybrid-ONT-Illumina      #Set the job name to "JobExample2"
#SBATCH --time=2-0:00:00               #Set the wall clock limit to 6hr and 30min
#SBATCH --ntasks=80
#SBATCH --ntasks-per-node=80          #Request 20 tasks/cores per node
#SBATCH --mem=2900G                     #Request 2GB per node
#SBATCH --partition=bigmem
#SBATCH --output=stdout.%j.txt          # save stdout to file
#SBATCH --error=stderr.%j.txt           # save stderr to file 

##SET OPTIONAL DIRECTIVES
#SBATCH --mail-type=END              #Send email on all job events
#SBATCH --mail-user=devonjboland@tamu.edu    #Send all emails to email_address

##############
# Created on: 03-16-21
# Created by: Devon J. Boland
##############

# Load Modules
ml GCC/9.3.0 SPAdes/3.14.1-Python-3.8.2

# Determine File Paths
workdir=/scratch/user/devonjboland/Canu-3TBNode-CanuAssembly
pe1=$workdir/IlluminaDataSet/Race_L_R1.fastq.gz
pe2=$workdir/IlluminaDataSet/Race_L_R2.fastq.gz
ONT=$workdir/ONTDataset/RaceL-barcode10-cell1-2.fastq.gz
CanuPolishedContigs=$workdir/Pilon-ContigPolishing/CANU-pilon/LRaceTest2_PILON.contigs.fasta
output1=$workdir/SPades/ONT+Illumina
output2=$workdir/SPades/ONT+Illumina+Canu-T
output3=$workdir/SPades/ONT+Illumina+Canu-U

# gzip ONT data
#gzip -c $workdir/ONTDataset/RaceL-barcode10-cell1-2.fastq > $workdir/ONTDataset/RaceL-barcode10-cell1-2.fastq.gz
# Has already been gzipped

# Run SPades Assembly on ONT and Illumina Data
spades.py -1 $pe1 -2 $pe2 --nanopore $ONT -output $output1 -t 80 -m 2900

# Run SPades Assembly on ONT and Illumina Data, with Canu contigs > Trusted
spades.py -1 $pe1 -2 $pe2 --nanopore $ONT --trusted-contigs $CanuPolishedContigs -output $output2 -t 80 -m 2900 

# Run SPades Assembly on ONT and Illumina Data, with Canu contigs > Untrusted
spades.py -1 $pe1 -2 $pe2 --nanopore $ONT --untrusted-contigs $CanuPolishedContigs -output $output3 -t 80 -m 2900 
