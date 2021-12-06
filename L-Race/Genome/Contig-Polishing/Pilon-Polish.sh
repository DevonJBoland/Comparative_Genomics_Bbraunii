#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=PilonPolish      #Set the job name to "JobExample2"
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

# Set Variables and Directories
workdir=/scratch/user/devonjboland/Canu-3TBNode-CanuAssembly
contigs=$workdir/Flye/wo-subassemblies/Flye.fasta
index=FlyeContigs
illuminareads=$workdir/Reads/IlluminaDataSet
pe1=$illuminareads/Race_L_R1.fastq
pe2=$illuminareads/Race_L_R2.fastq
output=Flye-Pilon.fasta
outputdir=$workdir/Flye

# Load Bowtie2 Module
ml GCC/7.3.0-2.30
ml OpenMPI/3.1.1
ml BWA/0.7.17
ml SAMtools/1.9
ml Pilon/1.23-Java-11

# Build HISAT2 Index Files of L race contigs
bwa index -p $index $contigs 

#Map Illumina paired end reads to L race contigs
bwa mem -t 40 $index $pe1 $pe2 | samtools sort -@ 40 -o Aligned.bam -
# Index Bam File
samtools index -b -@ 80 Aligned.bam Aligned.bam.bai

# Correct ONT Contigs with Pilon Using Illumina Mapped Reads
java -Xmx2900G -jar $EBROOTPILON/pilon.jar --threads 80 --genome $contigs --frags Aligned.bam --output $output --outdir $outputdir
