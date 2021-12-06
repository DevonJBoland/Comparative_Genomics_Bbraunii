#BSUB -L /bin/bash
#BSUB -J LongQCAnalysisTest1
#BSUB -o LongQCAnalysisTest1.%J
#BSUB -W 24:00
#BSUB -n 20
#BSUB -M 2700
#BSUB -R 'rusage[mem=2700]'
#BSUB -R 'span[ptile=20]'
#BSUB -e LongQCAnalysisTest1.stderr.txt

#BSUB -u devonjboland@tamu.edu
#BSUB -B -N


###############
# Created on: 01-08-2021
# Created by: Devon J. Boland
###############

module load LongQC/1.2.0-foss-2020a-Python-3.8.2
longqc=$EBROOTLONGQC/longQC.py

workdir=/scratch/user/devonjboland/CSH_ONT_ABLRaces_Bbraunii
cell1=$workdir/cell1
cell2=$workdir/cell2

# Obtain LongQC Analysis Metrics
for i in {08..10};
do
	python $longqc sampleqc -p 20 -x ont-rapid -o $workdir/LongQC-Analysis/cell1  $cell1/barcode$i/barcode$i.fastq
wait
done

for i in {08..10};
do
	python $longqc sampleqc -p 20 -x ont-rapid -o $workdir/LongQC-Analysis/cell2 $cell2/barcode$i/barcode$i.fastq
wait
done