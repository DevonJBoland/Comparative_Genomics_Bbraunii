#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=TE-Identification-Test    #Set the job name to "JobExample2"
#SBATCH --time=02-00:00:00               #Set the wall clock limit to 6hr and 30min
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=80          #Request 20 tasks/cores per node
#SBATCH --mem=2900G                     #Request 2GB per node
#SBATCH --partition=bigmem
#SBATCH --output=stdout.%j.txt          # save stdout to file
#SBATCH --error=stderr.%j.txt           # save stderr to file 

##SET OPTIONAL DIRECTIVES
#SBATCH --mail-type=ALL              #Send email on all job events
#SBATCH --mail-user=devonjboland@tamu.edu    #Send all emails to

# Variables
repeats='consensi.fa' # Repeat modeler generated library
algae='Bbraunii-Creinhardtii-transcripts.fa' # Reference algae to align too to remove transcript hits
genome='Arace-Nuclear-Genome-v1.fasta' # Genome of organism
Extracted_Sequences='Extracted_Sequences' # path to output directory for sequence extraction
extract_TE='/scratch/user/devonjboland/venv/TE_Model_Identification.py' # path to custom python script
venv='/scratch/user/devonjboland/venv'
weblogovenv='/scratch/user/devonjboland/weblogo' # Assign path to weblogo venv
weblogo='/scratch/user/devonjboland/weblogo/bin/weblogo' # Assign path to CL weblogo
filenames='Filtered-consensi_split_files'
threads=$SLURM_CPUS_PER_TASK

##########################################################################################
# Load BLAST
#ml GCC/9.3.0  OpenMPI/4.0.3  BLAST+/2.10.1 seqtk/1.3 SAMtools/1.10
# Blast initial RepeatLibrary against Algae Transcripts
# This was already done: cat Bbraunii-transcript.fa, Creinhardtii-transcript.fa > Algae-transcripts.fa
#makeblastdb -in ${algae} -dbtype nucl
#blastn -query ${repeats} -db ${algae} -outfmt 6 -out repeats_toremove.out -max_target_seqs 5 -evalue 1e-5 -num_threads ${threads}
# Remove any repeat models matching Bbraunii Brace and Creinhardtii Transcripts
#samtools faidx ${repeats}
#awk '{print $1}' repeats_toremove.out | sort | uniq > transcript-hits.txt
#awk '{print $1}' ${repeats}.fai | grep -v -f transcript-hits.txt > remove.txt
#seqtk subseq ${repeats} remove.txt > Filtered-${repeats}
#wait
##########################################################################################
##########################################################################################
# After this step, the resulting Filtered-${repeats} file will need to be broken into separate individual fasta files
# This was done on my personal computer using the following program, as it could not be
# installed on the super computer:
# https://pypi.org/project/split-fasta/
# Installed: pip install split-pasta
# Command: splitfasta filename.fasta
# Afterwards, the entire directory containing each individual fasta file was uploaded back to the cluster
##########################################################################################
##########################################################################################
# Identify up to 50 hits for each Repeat in FilteredRepeats.lib
# Create list containing names of all fasta files
#ls ${filenames} > fasta_list.txt
# Blast FilteredRepeats.lib against Genenome Assembly
#makeblastdb -in ${genome} -dbtype nucl
#mkdir BlastnResults
#cat fasta_list.txt | while read line
#do
#blastn -query ${filenames}/${line} -db ${genome} -outfmt 6 -evalue 1e-5 -max_hsps 1 -num_alignments 50 -num_threads ${threads} -out BlastnResults/${line}-blastn-results.tsv
#done

# Extract aligned sequences and create multi sequence fasta files per putative TE
#ml GCCcore/10.2.0 Python/3.8.6
#source ${venv}/bin/activate
#mkdir Extracted_Sequences
#cat fasta_list.txt | while read line
#do
#python3 ${extract_TE} ${genome} BlastnResults/${line}-blastn-results.tsv 125 Extracted_Sequences/${line}-extracted.fasta
#done
#deactivate

# Combine putative TE seqeunce with blastn hit extracted sequences
#ls Extracted_Sequences | sed 's/-extracted.fasta//g' > only_hits.txt # Create a list of the original putative TE fasta files that generated a hit to the genome
#cat only_hits.txt | while read line
#do
#cat ${filenames}/${line} Extracted_Sequences/${line}-extracted.fasta > Extracted_Sequences/${line}
#rm Extracted_Sequences/${line}-extracted.fasta
#done

# Load MAFFT and Align Extracted Fasta Sequences
#ml iccifort/2020.1.217  impi/2019.7.217  MAFFT/7.471-with-extensions
ml GCC/8.3.0  MUSCLE/3.8.1551
mkdir MSAFiles # create output directory
cat only_hits.txt | while read line # loop through all files
do
#mafft --thread $SLURM_CPUS_PER_TASK Extracted_Sequences/${line} > MSAFiles/${line}-msa.fasta
muscle -in Extracted_Sequences/${line} -out MSAFiles/${line}-msa.fasta
done
ml purge # purge all modules

# Load Python3 Module
ml GCCcore/10.2.0 Python/3.8.6
source ${weblogovenv}/bin/activate


# Loop through MSA files
mkdir WebLogoFiles
cat only_hits.txt | while read line # loop through all files
do
${weblogo} -f MSAFiles/${line}-msa.fasta -D fasta -o WebLogoFiles/${line}-weblogo.png -F png -A dna -t ${line}
done
deactivate
ml purge # purge all modules
