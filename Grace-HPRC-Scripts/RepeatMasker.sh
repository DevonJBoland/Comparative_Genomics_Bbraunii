#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=RepeatMasker-Lrace    #Set the job name to "JobExample2"
#SBATCH --time=01-00:00:00               #Set the wall clock limit to 6hr and 30min
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=80          #Request 20 tasks/cores per node
#SBATCH --mem=2900G                     #Request 2GB per node
#SBATCH --partition=bigmem
#SBATCH --output=stdout.%j.txt          # save stdout to file
#SBATCH --error=stderr.%j.txt           # save stderr to file 

##SET OPTIONAL DIRECTIVES
#SBATCH --mail-type=ALL              #Send email on all job events
#SBATCH --mail-user=devonjboland@tamu.edu    #Send all emails to

## Assign Variables ##
genome='Lrace-Nuclear-Genome-v1.fasta'
brace='Bbraunii_502_v2.1.transcript.fa'
outputdir='/scratch/user/devonjboland/Annotation/Lrace/Masking-NoTECuration'
models='consensi.fa'

## Run RepeatModeler via Singularity ##
# Format genome file as blastdb
singularity exec /scratch/data/bio/containers/tetools_1.4.sif BuildDatabase -engine ncbi -name ${mydb} ${genome}
# de novo compute repeat library for genome using RepeatModeler
#singularity exec /scratch/data/bio/containers/tetools_1.4.sif RepeatModeler -database ${mydb} -engine ncbi -pa $SLURM_CPUS_PER_TASK

## Load RepeatMasker Module & BLAST ##
ml GCC/9.3.0  OpenMPI/4.0.3  BLAST+/2.10.1 seqtk/1.3 SAMtools/1.10
# Remove Brace Transcript Hits from RepeatrtModeler Library #
makeblastdb -in ${brace} -dbtype nucl
blastn -task megablast -query ${models} -db ${brace} -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' -max_target_seqs 25 -culling_limit 2 -num_threads $SLURM_CPUS_PER_TASK -evalue 1e-10 -out transcript-hits.tsv
samtools faidx ${models}
awk '{print $1}' transcript-hits.tsv | sort | uniq > hits.txt
awk '{print $1}' ${models}.fai | grep -v -f hits.txt > RemovefromModels.txt 
seqtk subseq ${models} RemovefromModels.txt > Filtered-${models}
ml purge
# Mask Genome Using Custom RepeatLibrary with RepeatMasker 
ml GCC/9.3.0  OpenMPI/4.0.3 RepeatMasker/4.1.2-p1-HMMER
rmblast='/scratch/user/devonjboland/Annotation/rmblast-2.11.0/bin'
RepeatMasker -e ncbi -rmblast_dir ${rmblast} -pa $SLURM_CPUS_PER_TASK -nolow -lib Filtered-${models} -dir ${outputdir} -xsmall ${genome}
