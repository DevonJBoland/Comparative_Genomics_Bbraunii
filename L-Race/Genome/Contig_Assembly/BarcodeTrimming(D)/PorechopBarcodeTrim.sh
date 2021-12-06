#BSUB -L /bin/bash
#BSUB -J BarcodeTrimmingTest2
#BSUB -o BarcodeTrimmingTest2.txt
#BSUB -W 24:00
#BSUB -n 20
#BSUB -M 2700
#BSUB -R 'rusage[mem=2700]'
#BSUB -R 'span[ptile=20]'
#BSUB -e BarcodeTrimmingTest2.stderr.txt

#BSUB -u devonjboland@tamu.edu
#BSUB -B -N


###############
# Created on: 01-13-2021
# Created by: Devon J. Boland
###############

# Load relevant modules and assign module vairable
module load Porechop/0.2.4-intel-2019b-Python-3.7.4
porechop=$EBROOTPORECHOP/bin/porechop

# Set file locations
workdir=/scratch/user/devonjboland/CSH_ONT_ABLRaces_Bbraunii
cell1=$workdir/cell1
cell2=$workdir/cell2
mkdir $workdir/cell1trimmed
mkdir $workdir/cell2trimmed
cell1trimmed=$workdir/cell1trimmed
cell2trimmed=$workdir/cell2trimmed

# Trim barcode sequence from ONT seuquence data
# For cell1
for i in {08..10};
do
	mkdir $cell1trimmed/barcode$i
	cd $cell1trimmed/barcode$i
	$porechop -i $cell1/barcode$i/ -o cell1_barcode"$i"trimmed.fasta -t 20
	wait
done

# For cell2
for i in {08..10};
do
	mkdir $cell2trimmed/barcode$i
	cd $cell2trimmed/barcode$i
	$porechop -i $cell2/barcode$i/ -o cell2_barcode"$i"trimmed.fasta -t 20
	wait
done
