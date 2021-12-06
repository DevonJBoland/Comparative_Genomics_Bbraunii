#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=Arace-SPades     #Set the job name to "JobExample2"
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

## Set Variables ##
workdir='/scratch/user/devonjboland/A-Race/'
contigs='/scratch/user/devonjboland/A-Race/assembly-pilonx2.fasta'
pe1reads='/scratch/user/devonjboland/A-Race/Race_A_R1.fastq.gz'
pe2reads='/scratch/user/devonjboland/A-Race/Race_A_R2.fastq.gz'
spadesoutdir='/scratch/user/devonjboland/A-Race/SPades/'
spadesscaffolds='/scratch/user/devonjboland/A-Race/SPades/scaffolds.fasta'
pilonout3='/scratch/user/devonjboland/A-Race/SPades/spades-pilonx1'
pilonout4='/scratch/user/devonjboland/A-Race/SPades/spades-pilonx2'
threads=80

#Step1. Hybrid Assembly Illumina Data with Canu-Polished Contigs
#Load Spades#
ml GCC/8.2.0-2.31.1 OpenMPI/3.1.3 SPAdes/3.13.1
# Assembly Reads/Contigs
spades.py -m 2900 -t $threads -1 $pe1reads -2 $pe2reads --trusted-contigs $contigs -o $spadesoutdir
# Purge Module #
ml purge

#Step 2. Polish SPades Assembly with Pilon using Illumina Reads
#Load BWA/SAMTools/Pilon#
ml GCC/7.3.0-2.30
ml OpenMPI/3.1.1
ml BWA/0.7.17
ml SAMtools/1.9
ml Pilon/1.23-Java-11
#Polish Flye Assembly with Illumian Reads#
## Build HISAT2 Index Files of contigs##
bwa index $spadesscaffolds
## Align Illumina Reads to Indexed Flye Assembly
bwa mem -t $threads $spadesscaffolds $pe1reads $pe2reads | samtools sort -o reads.sorted3.bam -T reads.tmp -
samtools index -@ $threads reads.sorted3.bam
##Correct ONT Contigs with Pilon Using Illumina Mapped Reads##
java -Xmx2900G -jar $EBROOTPILON/pilon.jar --threads $threads --genome $spadesscaffolds --frags reads.sorted3.bam --output $pilonout3

#Step 3. Second Pilon Polishing Step
## Build HISAT2 Index Files of contigs##
bwa index ${pilonout3}.fasta
## Align Illumina Reads to Indexed Flye Assembly
bwa mem -t $threads ${pilonout3}.fasta $pe1reads $pe2reads | samtools sort -o reads.sorted4.bam -T reads.tmp -
samtools index -@ $threads reads.sorted4.bam
##Correct ONT Contigs with Pilon Using Illumina Mapped Reads##
java -Xmx2900G -jar $EBROOTPILON/pilon.jar --threads $threads --genome ${pilonout3}.fasta --frags reads.sorted4.bam --output $pilonout4
#Purge Modules#
ml purge

#Step 4. Quast Analysis of Scaffolds
## Load Modules ##
ml GCC/9.3.0
ml OpenMPI/4.0.3
ml QUAST/5.0.2-Python-3.8.2
quast=/sw/eb/sw/QUAST/5.0.2-foss-2020a-Python-3.8.2/bin/quast
#Quast Analysis of 
mkdir ${spadesoutdir}Quast
python $quast -t $threads -o ${spadesoutdir}Quast $spadesscaffolds ${pilonout3}.fasta ${pilonout4}.fasta
#Purge Modules
ml purge

#Step5. BUSCO Analysis of Assembly
#Load BUSCO Module#
module load GCC/8.3.0 OpenMPI/3.1.4 BUSCO/4.1.2-Python-3.7.4
##Assembly.fasta BUSCO Analysis##
genome_file=$spadesscaffolds
busco_mode='genome'                 # genome, transcriptome, proteins
busco_lineage='/scratch/user/devonjboland/odb10/lineages/chlorophyta_odb10'
lineage_name=$(basename $busco_lineage)
augustus_species='chlamydomonas'
augustus_model='partial'            # partial (default), intronless, complete, atleastone, exactlyone
out_prefix="out_busco4_${lineage_name}_${augustus_species}.SPades-Scaffolds"
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
busco --offline --in $genome_file --mode $busco_mode --augustus_species $augustus_species --cpu $threads --lineage $busco_lineage --out $out_prefix --augustus_parameters="--genemodel=$augustus_model"
## Pilonx3.fasta BUSCO Analysis##
genome_file=${pilonout3}.fasta
busco_mode='genome'                 # genome, transcriptome, proteins
busco_lineage='/scratch/user/devonjboland/odb10/lineages/chlorophyta_odb10'
lineage_name=$(basename $busco_lineage)
augustus_species='chlamydomonas'
augustus_model='partial'            # partial (default), intronless, complete, atleastone, exactlyone
out_prefix="out_busco4_${lineage_name}_${augustus_species}.SPades-Pilonx1"
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
busco --offline --in $genome_file --mode $busco_mode --augustus_species $augustus_species --cpu $threads --lineage $busco_lineage --out $out_prefix --augustus_parameters="--genemodel=$augustus_model"
## Pilonx4.fasta BUSCO Analysis##
genome_file=${pilonout4}.fasta
busco_mode='genome'                 # genome, transcriptome, proteins
busco_lineage='/scratch/user/devonjboland/odb10/lineages/chlorophyta_odb10'
lineage_name=$(basename $busco_lineage)
augustus_species='chlamydomonas'
augustus_model='partial'            # partial (default), intronless, complete, atleastone, exactlyone
out_prefix="out_busco4_${lineage_name}_${augustus_species}.SPades-Pilonx2"
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
busco --offline --in $genome_file --mode $busco_mode --augustus_species $augustus_species --cpu $threads --lineage $busco_lineage --out $out_prefix --augustus_parameters="--genemodel=$augustus_model"
