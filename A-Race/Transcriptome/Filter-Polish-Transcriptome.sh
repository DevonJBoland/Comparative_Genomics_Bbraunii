#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=Filter-Polish-Transcriptome-Arace   #Set the job name to "JobExample2"
#SBATCH --time=02-00:00:00               #Set the wall clock limit to 6hr and 30min
#SBATCH --ntasks=80
#SBATCH --ntasks-per-node=80          #Request 20 tasks/cores per node
#SBATCH --mem=2900G                     #Request 2GB per node
#SBATCH --partition=bigmem
#SBATCH --output=stdout.%j.txt          # save stdout to file
#SBATCH --error=stderr.%j.txt           # save stderr to file 

##SET OPTIONAL DIRECTIVES
#SBATCH --mail-type=ALL              #Send email on all job events
#SBATCH --mail-user=devonjboland@tapestatmu.edu    #Send all emails to


## Assign Variables ##
isoform=/sw/eb/sw/Trinity/2.12.0-foss-2020a-Python-3.8.2/trinityrnaseq-v2.12.0/util/align_and_estimate_abundance.pl
filter=/sw/eb/sw/Trinity/2.12.0-foss-2020a-Python-3.8.2/trinityrnaseq-v2.12.0/util/filter_low_expr_transcripts.pl
transcriptome='/scratch/user/devonjboland/A-Race/Transcriptome-Assembly/trinity-out/Trinity-GG.fasta'
matrix='/scratch/user/devonjboland/A-Race/Transcriptome-Assembly/trinity-out/rsem_outdir/RSEM.isoforms.results'
filtertranscriptome='/scratch/user/devonjboland/A-Race/Transcriptome-Assembly/trinity-out/Arace_High_Exp_Iso_Transcriptome_v1.fasta'
pe1='/scratch/user/devonjboland/A-Race/Genome-Assembly/Race_A_R1.fastq.gz'
pe2='/scratch/user/devonjboland/A-Race/Genome-Assembly/Race_A_R2.fastq.gz'
threads=80
memoryGB=2900
memoryMB=2900000

## Load Modules ##
#ml GCC/9.3.0  OpenMPI/4.0.3  Trinity/2.12.0-Python-3.8.2
#ml Bowtie2/2.4.1  SAMtools/1.10  RSEM/1.3.3
## Detect Low Expression Isoforms ##
#perl ${isoform} --thread_count ${threads} --transcripts ${transcriptome} --seqType fa --left ${pe1} --right ${pe2} --est_method RSEM --aln_method bowtie2 --trinity_mode --prep_reference --output_dir rsem_outdir
## Remove/Filter Low Expression Isoforms ##
#perl ${filter} --matrix ${matrix} --transcripts ${transcriptome} --highest_iso_only --trinity_mode
#ml purge

## Load Modules ##
ml  GCC/7.3.0-2.30  OpenMPI/3.1.1  CD-HIT/4.8.1
## Remove Redudant Transcripts ##
cd-hit-est -i ${filter-transcriptome} -o Arace_HighExp_Iso_EST_transcriptome_v1.fasta -c 0.95 -n 10 -M ${memoryMB} -T ${threads}
cd-hit-est -i ${transcriptome} -o Arace_EST_transcriptome_v1.fasta -c 0.95 -n 10 -M ${memoryMB} -T ${threads}
#ml purge

