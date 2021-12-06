#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=PacBio-SMRT       #Set the job name to "JobExample2"
#SBATCH --time=2-0:00:00               #Set the wall clock limit to 6hr and 30min
#SBATCH --ntasks=80
#SBATCH --ntasks-per-node=80          #Request 20 tasks/cores per node
#SBATCH --mem=2900G                     #Request 2GB per node
#SBATCH --partition=bigmem
#SBATCH --output=stdout.%j.txt          # save stdout to file
#SBATCH --error=stderr.%j.txt           # save stderr to file 

##SET OPTIONAL DIRECTIVES
#SBATCH --mail-type=ALL              #Send email on all job events
#SBATCH --mail-user=devonjboland@tamu.edu    #Send all emails to

workdir='/scratch/user/devonjboland/Brace-pacbio/'
bam2fastq='/sw/eb/sw/SMRT-Link/7.0.1.66975-cli-tools-only/smrtcmds/bin/bam2fastq'
bax2bam='/sw/eb/sw/SMRT-Link/7.0.1.66975-cli-tools-only/smrtcmds/bin/bax2bam'
AYHHT='/scratch/user/devonjboland/Brace-pacbio/AYHHT/'
h5='/scratch/user/devonjboland/Brace-pacbio/Pacbio_h5/'
threads=80

# Load SMRT-Link Tools Modul #
ml SMRT-Link/7.0.1.66975-cli-tools-only

# Extract all Reads from Compressed Archives #
#cd $AYHHT
#cat list.txt | while read line
#do
#	tar -xzf $line -C $AYHHT/uncompressed/
#done
 
# Convert .bax.h5 Files To .bam Files #
## AYHHT ##
cd $AYHHT/uncompressed/
$bax2bam -o AYHHTbam -f list-tgz.txt
## h5 ##
cd $h5
$bax2bam -o h5bam -f list.txt
## Had issues with movie name, I may have to use another command tool to extract movie names ##
# https://www.pacb.com/wp-content/uploads/SMRT_Tools_Reference_Guide_v700.pdf #

# Convert Bam Files to Fastq Files #
## AYHHT ##
cd $AYHHT/uncompressed/
cat bam.txt | while read line
do
	$bam2fastq -o AYHHT-reads $line
done

cat bam.txt | while read line
do
	$bam2fastq -o h5-reads $line
done