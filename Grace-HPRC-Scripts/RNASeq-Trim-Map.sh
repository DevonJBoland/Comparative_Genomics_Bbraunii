#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=RNASeq-Trim-Map    #Set the job name to "JobExample2"
#SBATCH --time=02-00:00:00               #Set the wall clock limit to 6hr and 30min
#SBATCH --ntasks=80
#SBATCH --ntasks-per-node=80          #Request 20 tasks/cores per node
#SBATCH --mem=2900G                     #Request 2GB per node
#SBATCH --partition=bigmem
#SBATCH --output=stdout.%j.txt          # save stdout to file
#SBATCH --error=stderr.%j.txt           # save stderr to file 

##SET OPTIONAL DIRECTIVES
#SBATCH --mail-type=ALL              #Send email on all job events
#SBATCH --mail-user=devonjboland@tamu.edu    #Send all emails to


## Set Variables ##
workdir='/scratch/user/devonjboland/L-Race/Transcriptome-Assembly'
rnaseq='/scratch/user/devonjboland/L-Race/Transcriptome-Assembly/RNA-Seq'
genome='/scratch/user/devonjboland/L-Race/Transcriptome-Assembly/Masked-Assembly/Lrace-Nuclear-Genome-v1-hard_masked.fasta'
memory=2900
threads=80

## Load Seqtk, Trimmomatic, FastQC Modules ##
#ml GCC/7.3.0-2.30  OpenMPI/3.1.1 seqtk/1.3
#module load Trimmomatic/0.39-Java-11
#ml FastQC/0.11.9-Java-11
## Deinterleave RNA seq reads ##
#cd ${rnaseq}
#ls *.fastq > list.txt
#cat list.txt | while read line
#do
#seqtk seq $line -1 | pigz -p ${threads} -c > ${line}_R1.fq.gz
#seqtk seq $line -2 | pigz -p ${threads} -c > ${line}_R2.fq.gz
## Trim RNA-Seq Data Using Trimmomatic ##
#java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads $threads ${line}_R1.fq.gz ${line}_R2.fq.gz -baseout ${line}_trim ILLUMINACLIP:TruSeq.fa:2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:100
## FastQC Analysis of Trimmed RNA seq PE Reads ##
#fastqc -o FastQC -t ${threads} ${line}_trim_1P ${line}_trim_2P ${line}_trim_1U ${line}_trim_2U
#done
#ml purge

## Load Hisat2 Module(s) ##
ml GCC/7.3.0-2.30  OpenMPI/3.1.1 SAMtools/1.9 HISAT2/2.2.0
## Index Bbraunii Genome ##
#cd ${workdir}
hisat2-build --threads ${threads} ${genome} ht2
## Map Trimmed RNA Seq PE Reads to Bbraunii Genome ##
cat ${rnaseq}/List.txt | while read line
do
hisat2 -p ${threads} -x ht2 -1 ${rnaseq}/${line}_trim_1P -2 ${rnaseq}/${line}_trim_2P -U ${rnaseq}/${line}_trim_1U,${rnaseq}/${line}_trim_2U | samtools sort -o ${line}_reads.sorted.bam -T reads.tmp -
done
samtools merge mapped-trim-reads.bam B.braunii_race.L_Day.00.fastq_reads.sorted.bam B.braunii_race.L_Day.03.fastq_reads.sorted.bam B.braunii_race.L_Day.07.fastq_reads.sorted.bam B.braunii_race.L_Day.14.fastq_reads.sorted.bam B.braunii_race.L_Day.21.fastq_reads.sorted.bam B.braunii_race.L_Day.28.fastq_reads.sorted.bam
samtools index -@ $threads mapped-trim-reads.bam
ml purge

## Load Trinity ##
ml GCC/8.3.0 OpenMPI/3.1.4 Trinity/2.10.0-Python-3.7.4
## Run Trinity in Genome Guided Mode ##
Trinity --genome_guided_bam  mapped-trim-reads.bam --max_memory ${memory}G --genome_guided_max_intron 10000 --CPU ${threads}
##### Done Running Trinity #####

if [ $* ]; then
    # check full-length reconstruction stats:

    ${TRINITY_HOME}/util/misc/illustrate_ref_comparison.pl __indiv_ex_sample_derived/refSeqs.fa trinity_out_dir/Trinity.fasta 90

    ./test_FL.sh --query trinity_out_dir/Trinity.fasta --target __indiv_ex_sample_derived/refSeqs.fa --no_reuse
fi
ml purge

## Assign Variables ##
isoform=/sw/eb/sw/Trinity/2.12.0-foss-2020a-Python-3.8.2/trinityrnaseq-v2.12.0/util/align_and_estimate_abundance.pl
filter=/sw/eb/sw/Trinity/2.12.0-foss-2020a-Python-3.8.2/trinityrnaseq-v2.12.0/util/filter_low_expr_transcripts.pl
transcriptome='/scratch/user/devonjboland/L-Race/Transcriptome-Assembly/Masked-Assembly/trinity_out_dir/Trinity-GG.fasta'
matrix='/scratch/user/devonjboland/L-Race/Transcriptome-Assembly/Masked-Assembly/trinity_out_dir/rsem_outdir/RSEM.isoforms.results'
filtertranscriptome='/scratch/user/devonjboland/L-Race/Transcriptome-Assembly/Masked-Assembly/trinity_out_dir/Lrace_High_Iso_Exp_Transcriptome_v2.fasta'
pe1='/scratch/user/devonjboland/L-Race/Genome-Assembly/Race_L_R1.fastq.gz'
pe2='/scratch/user/devonjboland/L-Race/Genome-Assembly/Race_L_R2.fastq.gz'
threads=80
memoryGB=2900
memoryMB=2900000

## Load Modules ##
ml GCC/9.3.0  OpenMPI/4.0.3  Trinity/2.12.0-Python-3.8.2
ml Bowtie2/2.4.1  SAMtools/1.10  RSEM/1.3.3
## Detect Low Expression Isoforms ##
perl ${isoform} --thread_count ${threads} --transcripts ${transcriptome} --seqType fq --left ${pe1} --right ${pe2} --est_method RSEM --aln_method bowtie2 --trinity_mode --prep_reference --output_dir rsem_outdir
## Remove/Filter Low Expression Isoforms ##
perl ${filter} --matrix ${matrix} --transcripts ${transcriptome} --highest_iso_only --trinity_mode
ml purge

## Load Modules ##
ml  GCC/7.3.0-2.30  OpenMPI/3.1.1  CD-HIT/4.8.1
## Remove Redudant Transcripts ##
cd-hit-est -i ${filtertranscriptome} -o Lrace_HighExp_Iso_EST_transcriptome_v2.fasta -c 0.95 -n 10 -M ${memoryMB} -T ${threads}
cd-hit-est -i ${transcriptome} -o Lrace_EST_transcriptome_v2.fasta -c 0.95 -n 10 -M ${memoryMB} -T ${threads}
