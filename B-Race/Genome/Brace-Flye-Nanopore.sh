#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=B-Race-Flye-Nanopore      #Set the job name to "JobExample2"
#SBATCH --time=2-0:00:00               #Set the wall clock limit to 6hr and 30min
#SBATCH --ntasks=80
#SBATCH --ntasks-per-node=80          #Request 20 tasks/cores per node
#SBATCH --mem=2900G                     #Request 2GB per node
#SBATCH --partition=bigmem
#SBATCH --output=stdout.%j.txt          # save stdout to file
#SBATCH --error=stderr.%j.txt           # save stderr to file 

##SET OPTIONAL DIRECTIVES
#SBATCH --mail-type=ALL              #Send email on all job events
#SBATCH --mail-user=devonjboland@tamu.edu    #Send all emails to

porechopout='/scratch/user/devonjboland/B-Race/Flye-Nanopore/trimmed-Brace-cell1-2-barcode09-ont.fastq.gz'
threads=80

#Step 1. Assemble ONT data using Flye
#Load Flye#
ml iccifort/2020.1.217  impi/2019.7.217
ml Flye/2.8.1-Python-3.8.2
#Assemble Trimmed ONT Reads#
cat min-overlap.txt | while read line
do
mkdir Flye-${line}
flye --nano-raw ${porechopout} -o Flye-${line} -t ${threads}
done
#Purge Modules#
ml purge
