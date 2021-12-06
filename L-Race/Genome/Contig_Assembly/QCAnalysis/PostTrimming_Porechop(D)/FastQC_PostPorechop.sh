#BSUB -L /bin/bash
#BSUB -J FastQC_PostPorechopTest1
#BSUB -o FastQC_PostPorechopTest1.txt
#BSUB -W 24:00
#BSUB -n 20
#BSUB -M 2700
#BSUB -R 'rusage[mem=2700]'
#BSUB -R 'span[ptile=20]'
#BSUB -e FastQC_PostPorechopTest1.stderr.txt

#BSUB -u devonjboland@tamu.edu
#BSUB -B -N


###############
# Created on: 01-15-2021
# Created by: Devon J. Boland
###############

# Load Modules and Assign Folders/Variables
module load FastQC/0.11.9-Java-1.8.0
fastqc=$EBROOTFASTQC/fastqc

workdir=/scratch/user/devonjboland/CSH_ONT_ABLRaces_Bbraunii
cell1trimmed=$workdir/cell1trimmed
cell2trimmed=$workdir/cell2trimmed
mkdir $workdir/PostTrim-FastQCAnalysis
outputdir=$workdir/PostTrim-FastQCAnalysis

# For cell1trimmed
for i in {08..10};
doß
	mkdir $outputdir/barcode$i
	$fastqc -o $outputdir/barcode$i -t 20 $cell1trimmed/barcode$i/cell1_barcode"$i"trimmed.fastq
wait
done

# For cell2trimmed
for i in {08..10};
do
	$fastqc -o $outputdir/barcode$i -t 20 $cell2trimmed/barcode$i/cell2_barcode"$i"trimmed.fastq
wait
done