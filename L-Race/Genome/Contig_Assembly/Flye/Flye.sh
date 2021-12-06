#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=FlyeAssemler        # job name
#SBATCH --time=02-00:00:00             # max job run time dd-hh:mm:ss
#SBATCH --ntasks=80
#SBATCH --ntasks-per-node=80          #Request 20 tasks/cores per node
#SBATCH --mem=2900G                    # total memory per node
#SBATCH --partition=bigmem          # use large (3TB) memory node
#SBATCH --output=stdout.%j.txt          # save stdout to file
#SBATCH --error=stderr.%j.txt           # save stderr to file

##############
# Created on: 03-30-21
# Created by: Devon J. Boland
##############

# Flye Assembler Test Script

# Load Modules
ml iccifort/2020.1.217  impi/2019.7.217 Flye/2.8.1-Python-3.8.2

# Vairables and Pathways
workdir=/scratch/user/devonjboland/Canu-3TBNode-CanuAssembly
ONTreads=$workdir/Reads/ONTDataset/RaceL-barcode10-cell1-2.fastq
CanuPolished=$workdir/Assemblies/CanuLrace-PILON.contigs.fasta
Outputdir1=$workdir/Flye/w-subassemblies
Outputdir2=$workdir/Flye/wo-subassemblies

# Execute Flye Assembler with Canu Contigs
flye --subassemblies $CanuPolished --out-dir $Outputdir1 --genome-size 211m -t 80

# Execute Flye Assembler without Canu Contigs
#flye --nano-raw $ONTreads --out-dir $Outputdir2 --genome-size 211m -t 80
