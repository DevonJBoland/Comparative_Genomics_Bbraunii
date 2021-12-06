#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=EDTA-Analysis    #Set the job name to "JobExample2"
#SBATCH --time=02-00:00:00               #Set the wall clock limit to 6hr and 30min
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=80          #Request 20 tasks/cores per node
#SBATCH --mem=2900G                     #Request 2GB per node
#SBATCH --partition=bigmem
#SBATCH --output=stdout.%j.txt          # save stdout to file
#SBATCH --error=stderr.%j.txt           # save stderr to file 

##SET OPTIONAL DIRECTIVES
#SBATCH --mail-type=ALL              #Send email on all job events
#SBATCH --mail-user=devonjboland@tamu.edu    #Send all emails to

# Assign Variables
workdir='/scratch/user/devonjboland/Repeat-Content-Analysis'
#aracedir=${workdir}/'Arace'
#araceg=${aracedir}/'Arace-Nuclear-Genome-v1-simple.fasta'
#aracecds=${aracedir}/'Arace-Predicted-Genes-CDS-v1.fasta'
#aracebed=${aracedir}/'sorted-Arace-gene.bed'
#bracedir=${workdir}/'Brace'
#braceg=${bracedir}/'Bbraunii_502_2.0-simple.fa'
#bracecds=${bracedir}/'Bbraunii_502_v2.1.cds.fa'
#bracebed=${bracedir}/'sorted.Brace.gene.bed'
#lracedir=${workdir}/'Lrace'
#lraceg=${lracedir}/'Lrace-Nuclear-Genome-v1-simple.fasta'
#lracecds=${lracedir}/'Lrace-Predicted-Genes-v1-CDS.fasta'
#lracebed=${lracedir}/'sorted.Lrace.gene.bed'
#chlredir=${workdir}/'Chlre'
#chlreg=${chlredir}/'Chlre5_6_AssemblyScaffolds_Repeatmasked.fasta'
#chlrecds=${chlredir}/'Chlre5_6_GeneCatalog_CDS_20200117.fasta'
#coccodir=${workdir}/'Cocco'
#coccog=${coccodir}/'Coccomyxa_C169_v2_masked_genomic_scaffolds.fasta'
cvarg="ChlNC64A_1_nuclear_scaffolds.fasta"
cvardir="/scratch/user/devonjboland/Annotation/Cvariabilis"
threads=$SLURM_CPUS_PER_TASK

# Use EDTA to Annotate and Predict TEs in Bbraunii
export SINGULARITY_BINDPATH="/scratch,$TMPDIR"
# Arace
#cd ${aracedir}
#singularity exec /sw/hprc/sw/bio/containers/edta_1.9.5.sif EDTA.pl --genome ${araceg} --cds ${aracecds} --overwrite 1 --anno 1 --evaluate 1 --threads ${threads}
# Brace
#cd ${bracedir}
#singularity exec /sw/hprc/sw/bio/containers/edta_1.9.5.sif EDTA.pl --genome ${braceg} --cds ${bracecds} --overwrite 1 --anno 1 --evaluate 1 --threads ${threads}
# Lrace
#cd ${lracedir}
#singularity exec /sw/hprc/sw/bio/containers/edta_1.9.5.sif EDTA.pl --genome ${lraceg} --cds ${lracecds} --overwrite 1 --anno 1 --evaluate 1 --threads ${threads}
# Chlamy
#cd ${chlredir}
#singularity exec /sw/hprc/sw/bio/containers/edta_1.9.5.sif EDTA.pl --genome ${chlreg} --cds ${chlrecds} --overwrite 1 --anno 1 --evaluate 1 --threads ${threads}
# Coccomyxia
#cd ${coccodir}
#singularity exec /sw/hprc/sw/bio/containers/edta_1.9.5.sif EDTA.pl --genome ${coccog} --overwrite 1 --anno 1 --evaluate 1 --threads ${threads}
cd ${cvardir}
singularity exec /sw/hprc/sw/bio/containers/edta_1.9.5.sif EDTA.pl --genome ${cvarg} --overwrite 1 --anno 1 --evaluate 1 --threads ${threads}
