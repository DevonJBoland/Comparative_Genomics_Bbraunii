#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=Nanopolish       #Set the job name to "JobExample2"
#SBATCH --time=2-0:00:00               #Set the wall clock limit to 6hr and 30min
#SBATCH --ntasks=80
#SBATCH --ntasks-per-node=80          #Request 20 tasks/cores per node
#SBATCH --mem=2900G                     #Request 2GB per node
#SBATCH --partition=bigmem
#SBATCH --output=stdout.%j.txt          # save stdout to file
#SBATCH --error=stderr.%j.txt           # save stderr to file 

##SET OPTIONAL DIRECTIVES
#SBATCH --mail-type=ALL              #Send email on all job events
#SBATCH --mail-user=devonjboland@tamu.edu    #Send all emails to email_address

# Set Variables and Directories
workdir=/scratch/user/devonjboland/Canu-3TBNode-CanuAssembly
assembly=$workdir/Nanopolish/Flye.fasta
ONTData=$workdir/Reads/ONTDataset
reads=$ONTData/RaceL-barcode10-cell1-2.fastq
fast5PAF=$ONTData/barcode10_fast5/PAF/
fast5PAG=$ONTData/barcode10_fast5/PAF/
summarytxt=$ONTData/barcode10_fast5/
scripts=/sw/eb/sw/nanopolish/0.13.1-foss-2018b-Python-3.6.6/scripts
logfile=log.txt

# Load Modules
ml GCC/7.3.0-2.30
ml OpenMPI/3.1.1
ml BWA/0.7.17
ml SAMtools/1.9
ml nanopolish/0.13.1-Python-3.6.6
ml parallel/20180822

#Index the fastq files
nanopolish index -d $fast5PAF -d $fast5PAG -f $summarytxt $reads
echo "Fast 5 Indexing DONE" >> $logfile
# Index the draft genome
bwa index -p index $assembly 
echo "bwa Indexing DONE" >> $logfile
# Align the basecalled reads to the draft sequence
bwa mem -x index -t 80 $assembly $reads | samtools sort -o reads.sorted.bam -T reads.tmp -
samtools index -@ 80 reads.sorted.bam
echo "Align reads to assembly DONE" >> $logfile
# Correct Canu Contigs (polish) with Nanopolish
python3 $scripts/nanopolish_makerange.py $assembly | parallel --results nanopolish.results -P 40 nanopolish variants --consensus -o polished.{1}.vcf -w {1} -r $reads -b reads.sorted.bam -g $assembly -t 40 --min-candidate-frequency 0.1
echo "Correct individual contigs DONE" >> $logfile
# Merge all blocks into a single fasta file
nanopolish vcf2fasta -g $assembly polished.*.vcf > Flye-Nanopolish.fasta
echo "Merge all blocks back into a single file DONE" >> $logfile
