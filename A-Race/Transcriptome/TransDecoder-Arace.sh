#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=TransDecoder-Arace    #Set the job name to "JobExample2"
#SBATCH --time=02-00:00:00               #Set the wall clock limit to 6hr and 30min
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=80          # CPUs (threads) per command
#SBATCH --mem=2900G                     #Request 2GB per node
#SBATCH --partition=bigmem
#SBATCH --output=stdout.%j.txt          # save stdout to file
#SBATCH --error=stderr.%j.txt           # save stderr to file 

##SET OPTIONAL DIRECTIVES
#SBATCH --mail-type=ALL              #Send email on all job events
#SBATCH --mail-user=devonjboland@tapestatmu.edu    #Send all emails to

## Set Variables ##
transcripts='Arace_EST_transcriptome_v1.fasta'
longorfs='/sw/eb/sw/TransDecoder/5.5.0-foss-2020b/TransDecoder.LongOrfs'
predict='/sw/eb/sw/TransDecoder/5.5.0-foss-2020b/TransDecoder.Predict'
nr="/scratch/data/bio/blastdb-2021.05.08/nr"
pfamA='/scratch/user/devonjboland/Cece-Research/Pfam-A.hmm'
threads=$SLURM_CPUS_PER_TASK

## Load Modules ##
#ml GCC/10.2.0  OpenMPI/4.0.5 TransDecoder/5.5.0
## Detect Long ORF ##
#${longorfs} -t ${transcripts}
#cp *transdecoder_dir/longest_orfs.pep ../transdecoder_outdir/

#ml purge

## Load Blast Modules ##
#ml  GCC/10.2.0  OpenMPI/4.0.5 BLAST+/2.11.0 
## Blastp LongORF Against nr Database ##
#blastp -query longest_orfs.pep -db ${nr} -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads ${threads} > blastp.outfmt6

#ml purge

## Load HMMER Modules ##
#ml GCC/10.2.0  OpenMPI/4.0.5  HMMER/3.3.2
## Perform HMMSEARCH Against PFAM-A Database ##
#hmmscan --cpu ${threads} --domtblout pfam.domtblout ${pfamA} longest_orfs.pep

#ml purge

## Load TransDecoder Modeules ##
ml GCC/10.2.0  OpenMPI/4.0.5 TransDecoder/5.5.0
## Predict ORFS ##
${predict}  -t ${transcripts} --retain_pfam_hits pfam.domtblout
