#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=Blobtools-Arace      #Set the job name to "JobExample2"
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
assembly='/scratch/user/devonjboland/blobtools/A-Race/Arace-Quickmerge-ml5000-Pilon.fasta'
pe1reads='/scratch/user/devonjboland/A-Race/Race_A_R1.fastq.gz'
pe2reads='/scratch/user/devonjboland/A-Race/Race_A_R2.fastq.gz'
db='/scratch/data/bio/blastdb-2021.05.08/nt'
blastout='Arace-blast-hits'
threads=80

# Map PE Reads to Assembly
#Load BWA/SAMTools/Pilon#
#ml GCC/8.3.0
#ml Bowtie2/2.3.5.1
#ml SAMtools/1.10
## Build BWA Index Files of L race contigs##
#bowtie2-build --threads $threads -f $assembly genome
## Align Illumina Reads to Indexed Flye Assembly
#bowtie2 -p $threads -x genome -1 $pe1reads -2 $pe2reads | samtools sort -o reads.sorted.bam -T reads.tmp -
#samtools index -@ $threads reads.sorted.bam
## Purge Modules ##
#ml purge

# Blast Assembly Against NCBI nt Databse
## Load Modules ##
ml GCC/10.2.0  OpenMPI/4.0.5 BLAST+/2.11.0 
## Perform blastn ##
blastn -task megablast -query $assembly -db ${db} -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' -max_target_seqs 1 -max_hsps 1 -num_threads $threads -evalue 1e-25 -out $blastout
## Purge Modules ##
ml purge

# Identify and Map Contaminant DNA via Blobtools #
## Load Modules ##
ml GCC/7.3.0-2.30  OpenMPI/3.1.1 BlobTools/20180528-Python-2.7.15
# Convert BAM sorted File to lower size Cov File #
blobtools map2cov -o reads.sorted.cov -i $assembly -b reads.sorted.bam  
# Run Blobtools to create blob plot, txt file, etc.
blobtools create -i ${assembly} -c reads.sorted.cov.reads.sorted.bam.cov -t ${blastout}-tax.${blastout}.taxified.out -o A-race-blob && \
blobtools view -r order -i A-race-blob.blobDB.json && \
blobtools plot -r order -i A-race-blob.blobDB.json
