#BSUB -L /bin/bash
#BSUB -J MasurcaAssembly
#BSUB -o stdout.txt
#BSUB -e stderr.txt
#BSUB -n 40
#BSUB -R "select[mem2tb]"
#BSUB -R "span[ptile=40]"
#BSUB -R "rusage[mem=49750]"
#BSUB -M 49750
#BSUB -q xlarge
#BSUB -W 240:00
	
#BSUB -u devonjboland@tamu.edu
#BSUB -B -N
	
#############################
# Created on: 03-16-2021
# Created by: Devon J. Boland
#############################
	
# Load Westmere module for Large >1TB memory jobs
module load Westmere
module load MaSuRCA/3.4.1-foss-2018b-Perl-5.28.0
	
# Execute Assemble.sh
	
./Assemble.sh
