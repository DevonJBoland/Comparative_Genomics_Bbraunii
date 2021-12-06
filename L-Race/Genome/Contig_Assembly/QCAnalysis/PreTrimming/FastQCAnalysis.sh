#BSUB -L /bin/bash
#BSUB -J FastQCAnalysisTest1
#BSUB -o FastQCAnalysisTest1.%J
#BSUB -W 24:00
#BSUB -n 20
#BSUB -M 2700
#BSUB -R 'rusage[mem=2700]'
#BSUB -R 'span[ptile=20]'
#BSUB -e FastQCAnalysisTest1.stderr.txt

#BSUB -u devonjboland@tamu.edu
#BSUB -B -N


###############
# Created on: 01-08-2021
# Created by: Devon J. Boland
###############

module load FastQC/0.11.9-Java-1.8.0
fastqc=$EBROOTFASTQC/fastqc

workdir=/scratch/user/devonjboland/CSH_ONT_ABLRaces_Bbraunii
cell1=$workdir/cell1
cell2=$workdir/cell2

# Collect Quality Metrics on ONT Data
mkdir FastQC-Analysis
mkdir FastQC-Analysis/cell1
mkdir FastQC-Analysis/cell2
cell1out=FastQC-Analysis/cell1
cell2out=FastQC-Analysis/cell2

for i in {08..10};
do
	cat $cell1/barcode$i/*.fastq >> $cell1/barcode$i/barcode$i.fastq
	$fastqc -o $cell1out -t 20 $cell1/barcode$i/barcode$i.fastq
done
wait

for i in {08..10};
do
	cat $cell2/barcode$i/*.fastq >> $cell2/barcode$i/barcode$i.fastq
	$fastqc -o $cell2out -t 20 $cell2/barcode$i/barcode$i.fastq
done
