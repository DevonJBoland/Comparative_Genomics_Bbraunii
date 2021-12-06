#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=ONTContigs-Pilon      #Set the job name to "JobExample2"
#SBATCH --time=2-0:00:00               #Set the wall clock limit to 6hr and 30min
#SBATCH --ntasks=80
#SBATCH --ntasks-per-node=80          #Request 20 tasks/cores per node
#SBATCH --mem=2900G                     #Request 2GB per node
#SBATCH --partition=bigmem
#SBATCH --output=ONTContigs-Pilon.stdout.%j.txt          # save stdout to file
#SBATCH --error=ONTContigs-Pilon.stderr.%j.txt           # save stderr to file

##SET OPTIONAL DIRECTIVES
##SBATCH --mail-type=ALL              #Send email on all job events
##SBATCH --mail-user=devonjboland@tamu.edu    #Send all emails to email_address

# Load Modules
ml GCC/7.3.0-2.30
ml OpenMPI/3.1.1
ml SAMtools/1.9
ml Pilon/1.23-Java-11

# Set Variables
workdir=/scratch/user/devonjboland/Canu-3TBNode-CanuAssembly
contigs=$workdir/Canu-Assembled-Contigs/LraceTest2.contigs.fasta
frags=$workdir/Pilon-ContigPolishing/HISAT2/IlluminaAligned-ONTContigs/IlluminaAligned-Sorted.bam
outputdir=$workdir/Canu-Assembled-Contigs
output=LraceTest2_POSLISHED.contigs

# Index Bam File
samtools index -b -@ 80 $frags $frags.bai

# Correct ONT Contigs with Pilon Using Illumina Mapped Reads
java -Xmx2900G -jar $EBROOTPILON/pilon.jar --threads 20 --genome $contigs --frags $frags --output $output --outdir $outputdir