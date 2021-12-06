#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=Blast-Lrace      #Set the job name to "JobExample2"
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

# Create Variables #
workdir='/scratch/user/devonjboland/'
pe1reads='/scratch/user/devonjboland/L-Race/Race_L_R1.fastq.gz'
pe2reads='/scratch/user/devonjboland/L-Race/Race_L_R2.fastq.gz'
assembly='/scratch/user/devonjboland/blobtools/L-Race/assembly-pilonx1.fasta'
db='/scratch/data/bio/blastdb-2021.05.08/nt'
blastout='Lrace-blast-hits'
threads=80

# Map PE Reads to Assembly
#Load BWA/SAMTools/Pilon#
#ml GCC/7.3.0-2.30
#ml OpenMPI/3.1.1
#ml BWA/0.7.17
#ml SAMtools/1.9
## Build BWA Index Files of L race contigs##
#bwa index -p assembly $assembly
## Align Illumina Reads to Indexed Flye Assembly
#bwa mem -t $threads assembly $pe1reads $pe2reads | samtools sort -o reads.sorted.bam -T reads.tmp -
#samtools index -@ $threads reads.sorted.bam
## Purge Modules ##
#ml purge

# Blast Assembly Against NCBI nt Databse
## Load Modules ##
ml GCC/9.3.0  OpenMPI/4.0.3 BLAST+/2.10.1 
## Perform blastn ##
blastn -task megablast -query ${assembly} -db ${db} -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' -max_target_seqs 1 -max_hsps 1 -num_threads $threads -evalue 1e-25 -out $blastout
## Purge Modules ##
ml purge

# Identify and Map Contaminant DNA via Blobtools #
## Load Modules ##
ml GCC/7.3.0-2.30  OpenMPI/3.1.1 BlobTools/20180528-Python-2.7.15
# Convert BAM sorted File to lower size Cov File #
#blobtools map2cov -o map2cov-bam.reads.sorted.bam.cov -i $assembly -b reads.sorted.bam
# Run Blobtools to create blob plot, txt file, etc.
blobtools create -i ${assembly} -c map2cov-bam.reads.sorted.bam.cov -t ${blastout} -o L-race-blob && \
blobtools view -r order -i L-race-blob.blobDB.json && \
blobtools plot -r order -i L-race-blob.blobDB.json
