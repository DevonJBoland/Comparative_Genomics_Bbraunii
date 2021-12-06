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


##############
# Created on: 03-26-21
# Created by: Devon J. Boland
##############

# Load Modules
ml GCC/7.3.0-2.30
ml OpenMPI/3.1.1
ml HISAT2/2.2.0
ml SAMtools/1.9
ml Pilon/1.23-Java-11

# Set Variable Paths
workdir=/scratch/user/devonjboland/Canu-3TBNode-CanuAssembly
assemblies=$workdir/Assemblies/Merged-Canu-Masurca-Lrace.ml500000.scaffolds.fasta
illuminareads=$workdir/Reads/IlluminaDataSet
pe1=$illuminareads/Race_L_R1.fastq
pe2=$illuminareads/Race_L_R2.fastq
ht2index=Index
frags=IlluminaAligned-Sorted.bam
output=Merged-Canu-Masurca-Lrace.ml500000.scaffolds-PILON.fasta
outputdir=$workdir/Assemblies

# Build HISAT2 Index Files of L race contigs
hisat2-build -p 80 $assemblies $ht2index

#Map Illumina paired end reads to L race contigs
hisat2  --sensitive -p 40 -x $ht2index -1 $pe1 -2 $pe2 | samtools view -@ 40 -bSh - -o IlluminaAligned.bam

# Sort BAM file and Index BAM file
samtools sort -o IlluminaAligned-Sorted.bam -O bam -@ 80 IlluminaAligned.bam
samtools index -b -@ 80 IlluminaAligned-Sorted.bam IlluminaAligned-Sorted.bam.bai

# Correct ONT Contigs with Pilon Using Illumina Mapped Reads
java -Xmx2900G -jar $EBROOTPILON/pilon.jar --threads 80 --genome $assemblies --frags $frags --output $output --outdir $outputdir