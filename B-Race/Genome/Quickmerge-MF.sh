#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=Brace-Quickmerge     #Set the job name to "JobExample2"
#SBATCH --time=5:00:00               #Set the wall clock limit to 6hr and 30min
#SBATCH --ntasks=80
#SBATCH --ntasks-per-node=80          #Request 20 tasks/cores per node
#SBATCH --mem=2900G                     #Request 2GB per node
#SBATCH --partition=bigmem
#SBATCH --output=stdout.%j.txt          # save stdout to file
#SBATCH --error=stderr.%j.txt           # save stderr to file 

##SET OPTIONAL DIRECTIVES
#SBATCH --mail-type=ALL              #Send email on all job events
#SBATCH --mail-user=devonjboland@tamu.edu    #Send all emails to

## Set Variables ##
pe1reads='/scratch/user/devonjboland/B-Race/SXPX_150_R1.fastq.gz'
pe2reads='/scratch/user/devonjboland/B-Race/SXPX_150_R2.fastq.gz'
threads=80


#Step 1. Polish Flye/Masurca Assemblies with Pilon using Illumina Reads
#Load BWA/SAMTools/Pilon#
#ml GCC/7.3.0-2.30
#ml OpenMPI/3.1.1
#ml BWA/0.7.17
#ml SAMtools/1.9
#ml Pilon/1.23-Java-11
#Polish Flye Assembly with Illumian Reads#
#ls *.fasta > assembly-list.txt
#cat assembly-list.txt | while read line
#do
## Build HISAT2 Index Files of contigs ##
#bwa index ${line}
## Align Illumina Reads to Indexed Flye Assembly
#bwa mem -t $threads ${line} ${pe1reads} ${pe2reads} | samtools sort -o ${line}-reads.sorted.bam -T reads.tmp -
#samtools index -@ $threads ${line}-reads.sorted.bam
## Correct ONT Contigs with Pilon Using Illumina Mapped Reads ##
#java -Xmx2900G -jar $EBROOTPILON/pilon.jar --threads $threads --genome ${line} --frags ${line}-reads.sorted.bam --output ${line}-pilon
#done
#ml purge

#Step 2. Merge Flye and Masurca Assemblies - Loop
# Load Modueles 
ml GCC/9.3.0 quickmerge/0.3
# Loop through ml flag 5000 - 500000 by power of 10 (total of three iterations)
cat ml-list.txt | while read line
do
merge_wrapper.py -pre brace-merged.ml${line} -l 173869 -ml ${line} brace-masurca.fasta-pilon.fasta brace-flye.fasta-pilon.fasta
done
ml purge

#Step 3. Second Pilon Polishing Step
## Load Modules ## 
ml GCC/7.3.0-2.30
ml OpenMPI/3.1.1
ml BWA/0.7.17
ml SAMtools/1.9
ml Pilon/1.23-Java-11
## Loop through merged assemblies and polish ##
ls brace-merged.ml* > quickmerge-list.txt
cat quickmerge-list.txt | while read line
do
## Build HISAT2 Index Files of contigs##
bwa index ${line}
## Align Illumina Reads to Indexed Flye Assembly
bwa mem -t $threads ${line} $pe1reads $pe2reads | samtools sort -o ${line}-reads.sorted.bam -T reads.tmp -
samtools index -@ $threads ${line}-reads.sorted.bam
##Correct ONT Contigs with Pilon Using Illumina Mapped Reads##
java -Xmx2900G -jar $EBROOTPILON/pilon.jar --threads $threads --genome ${line} --frags ${line}-reads.sorted.bam --output ${line}.pilon
done
ml purge

#Step 4. Quast Analysis of Scaffolds
## Load Modules ##
ml GCC/9.3.0
ml OpenMPI/4.0.3
ml QUAST/5.0.2-Python-3.8.2
quast=/sw/eb/sw/QUAST/5.0.2-foss-2020a-Python-3.8.2/bin/quast
#Quast Analysis of 
mkdir Quast
python $quast -t $threads -o Quast *.fasta
