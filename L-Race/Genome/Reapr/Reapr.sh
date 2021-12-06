#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=Reapr      #Set the job name to "JobExample2"
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


# Load Modules
ml GCC/8.2.0-2.31.1  OpenMPI/3.1.3 Reapr/1.0.18

# Assembly
workdir=/scratch/user/devonjboland/Canu-3TBNode-CanuAssembly
assembly=$workdir/Pilon-ContigPolishing/Merged-pilon/Merged+Canu_PILONx2.scaffolds.fasta
new_assembly=LraceAssembly-Reapr
# Reads
illuminareads=$workdir/IlluminaDataSet
pe1=$illuminareads/Race_L_R1.fastq
pe2=$illuminareads/Race_L_R2.fastq
# Output Directory
outputdir=/scratch/user/devonjboland/Canu-3TBNode-CanuAssembly/Reapr

# Run facheck on assembly and rename to new assembly
reapr facheck $assembly $new_assembly

# Align Illumina data and run pipeline (4.4 in Reapr manual)
reapr perfectmap $new_assembly.fa $pe1 $pe2 270 perfect
reapr smaltmap -n 80 $new_assembly.fa $pe1 $pe2 mapped.bam
reapr pipeline $new_assembly.fa mapped.bam $outputdir perfect
