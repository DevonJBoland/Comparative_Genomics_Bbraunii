#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=B-Race-Flye-SubAssembly      #Set the job name to "JobExample2"
#SBATCH --time=2-0:00:00               #Set the wall clock limit to 6hr and 30min
#SBATCH --ntasks=80
#SBATCH --ntasks-per-node=80          #Request 20 tasks/cores per node
#SBATCH --mem=2900G                     #Request 2GB per node
#SBATCH --partition=bigmem
#SBATCH --output=B-race-stdout.%j.txt          # save stdout to file
#SBATCH --error=B-race-stderr.%j.txt           # save stderr to file 

##SET OPTIONAL DIRECTIVES
#SBATCH --mail-type=ALL              #Send email on all job events
#SBATCH --mail-user=devonjboland@tamu.edu    #Send all emails to

## Set Variables##
quickmerge='/scratch/user/devonjboland/B-Race/Flye-Sub/Brace-Quickmerge-ml5000-Pilon.fasta'
bracev1='/scratch/user/devonjboland/B-Race/Flye-Sub/Bbraunii_502_2.0.fa'
flye='/scratch/user/devonjboland/B-Race/Flye-Sub/brace-flye.fasta-pilon.fasta'
pe1reads='/scratch/user/devonjboland/B-Race/Brace_SXPX_R1.fastq.gz'
pe2reads='/scratch/user/devonjboland/B-Race/Brace_SXPX_R2.fastq.gz'
pilonout='Brace-FlyeSub-Pilon'
threads=80

#Step 1. Assemble ONT data using Flye
#Load Flye#
ml iccifort/2020.1.217  impi/2019.7.217
ml Flye/2.8.1-Python-3.8.2
#Assemble Trimmed ONT Reads#
cat min-overlap.txt | while read line
do
mkdir Flye-${line}
flye --subassemblies ${quickmerge} ${bracev1} ${flye} -t ${threads} -i 0 -m ${line} -o Flye-${line}
#Purge Modules#
ml purge

# Step 2. Map PE Reads to Assembly #
#Load BWA/SAMTools/Pilon#
ml GCC/8.3.0
ml Bowtie2/2.3.5.1
ml SAMtools/1.10
## Build Bowtie2Index Files of L race contigs##
bowtie2-build --threads $threads -f assembly.fasta genome
## Align Illumina Reads to Indexed Flye Assembly
bowtie2 -p $threads -x genome -1 $pe1reads -2 $pe2reads | samtools sort -o reads.sorted.bam -T reads.tmp -
samtools index -@ $threads reads.sorted.bam
## Purge Modules ##
ml purge

# Step 3. Polish Assembly.fasta #
# Load Modules #
ml GCC/7.3.0-2.30
ml Pilon/1.23-Java-11
# Pilon Polish #
java -Xmx2900G -jar $EBROOTPILON/pilon.jar --threads ${threads} --genome assembly.fasta --frags reads.sorted.bam --output ${pilonout}
# Purge Modules #
ml purge

# Step 4. Quast Analysis #
# Load Modules #
ml GCC/9.3.0
ml OpenMPI/4.0.3
ml QUAST/5.0.2-Python-3.8.2
quast=/sw/eb/sw/QUAST/5.0.2-foss-2020a-Python-3.8.2/bin/quast
mkdir Quast
python $quast -t $threads -o Quast ${pilonout}.fasta
# Purge Modules #
ml purge

# Step 5. BUSCO Analysis #
# Load Modules #
module load GCC/10.2.0 OpenMPI/4.0.5 BUSCO/5.1.2
genome_file=${pilonout}.fasta
######## PARAMETERS ########
threads=$SLURM_CPUS_PER_TASK
busco_mode='genome'                 # genome, transcriptome, proteins
busco_lineage='/scratch/user/devonjboland/odb10/lineages/chlorophyta_odb10'         # check available databases at /scratch/data/bio/busco5/lineages/
augustus_species='chlamydomonas' # see available species here: /sw/eb/sw/AUGUSTUS/3.4.0-foss-2020b/config/species/
augustus_model='partial'            # partial (default), intronless, complete, atleastone, exactlyone
########## OUTPUTS #########
out_prefix="out_busco5_chlorophyta_odb10_chlamydomonas.Brace-FlyeSub-Pilon"
########## CONFIG ##########
# this section will copy augustus config to $SCRATCH;  no need to change this section
if [ ! -d "$SCRATCH/my_augustus_config/config" ]; then
  echo "Copying AUGUSTUS config directories to $SCRATCH/my_augustus_config"
  mkdir $SCRATCH/my_augustus_config
    if [ "-$EBROOTAUGUSTUS" == "-" ]; then
        echo "Augustus module not loaded"; exit 1
    fi
  rsync -r $EBROOTAUGUSTUS/ $SCRATCH/my_augustus_config
  chmod -R 755 $SCRATCH/my_augustus_config
fi
export AUGUSTUS_CONFIG_PATH="$SCRATCH/my_augustus_config/config"
################################### COMMANDS ###################################
busco --offline --in $genome_file --mode $busco_mode --cpu $threads --lineage $busco_lineage --out $out_prefix --augustus --augustus_species $augustus_species --augustus_parameters="--genemodel=$augustus_model"
################################################################################
