#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=B-Race-Assembly      #Set the job name to "JobExample2"
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
workdir='/scratch/user/devonjboland/B-Race'
pe1reads='/scratch/user/devonjboland/B-Race/Brace_SXPX_R1.fastq.gz'
pe2reads='/scratch/user/devonjboland/B-Race/Brace_SXPX_R2.fastq.gz'
ontreads='/scratch/user/devonjboland/B-Race/Brace-cell1-2-barcode09-ont.fastq.gz'
porechopout='/scratch/user/devonjboland/B-Race/trimmed-Brace-cell1-2-barcode09-ont.fastq.gz'
pilonout1='/scratch/user/devonjboland/B-Race/assembly-pilonx1'
pilonout2='/scratch/user/devonjboland/B-Race/assembly-pilonx2'
threads=80

pigz -c -p 80 Brace_SXPX_1P > Brace_SXPX_1P.fastq.gz
pigz -c -p 80 Brace_SXPX_2P > Brace_SXPX_2P.fastq.gz

#Step 1. Trim Barcodes/Adapters from ONT Reads
#Load Porechop#
#ml iccifort/2019.5.281  impi/2018.5.288
#ml Porechop/0.2.4-Python-3.7.4
#porechop=$EBROOTPORECHOP/bin/porechop
#Trim ONT Reads#
#$porechop -i $ontreads -o $porechopout -t $threads --end_threshold 90
#Purge Modules#
#ml purge

#Step 2. Assemble ONT data using Flye
#Load Flye#
#ml iccifort/2020.1.217  impi/2019.7.217
#ml Flye/2.8.1-Python-3.8.2
#Assemble Trimmed ONT Reads#
#flye --nano-raw $porechopout -o $workdir -t $threads
#Purge Modules#
#ml purge

#Step 3. Polish Flye Assembly Twice with Illumina Reads 
#Load BWA/SAMTools/Pilon#
ml GCC/7.3.0-2.30
ml OpenMPI/3.1.1
ml BWA/0.7.17
ml SAMtools/1.9
ml Pilon/1.23-Java-11
#Polish Flye Assembly with Illumian Reads#
## Build HISAT2 Index Files of L race contigs##
bwa index assembly.fasta
## Align Illumina Reads to Indexed Flye Assembly
bwa mem -t $threads assembly.fasta $pe1reads $pe2reads | samtools sort -o reads.sorted1.bam -T reads.tmp -
samtools index -@ $threads reads.sorted1.bam
##Correct ONT Contigs with Pilon Using Illumina Mapped Reads##
java -Xmx2900G -jar $EBROOTPILON/pilon.jar --threads $threads --genome assembly.fasta --frags reads.sorted1.bam --output $pilonout1

#Step 4. Second Pilon Polishing Step
## Build HISAT2 Index Files of L race contigs##
bwa index ${pilonout1}.fasta
## Align Illumina Reads to Indexed Flye Assembly
bwa mem -t $threads ${pilonout1}.fasta $pe1reads $pe2reads | samtools sort -o reads.sorted2.bam -T reads.tmp -
samtools index -@ $threads reads.sorted2.bam
##Correct ONT Contigs with Pilon Using Illumina Mapped Reads##
java -Xmx2900G -jar $EBROOTPILON/pilon.jar --threads $threads --genome ${pilonout1}.fasta --frags reads.sorted2.bam --output $pilonout2
#Purge Modules#
ml purge

#Step 5. Quast Analysis of Assembly
#Load Quast Module#
#ml GCC/9.3.0
#ml OpenMPI/4.0.3
#ml QUAST/5.0.2-Python-3.8.2
#quast=/sw/eb/sw/QUAST/5.0.2-foss-2020a-Python-3.8.2/bin/quast
#Quast Analysis of Flye Assembly, Pilonx1 Assembly, Pilonx2 Assembly
#mkdir $workdir/Quast
#python $quast -t $threads -o $workdir/Quast/ assembly.fasta ${pilonout1}.fasta ${pilonout2}.fasta
#Purge Modules
#ml purge

#Step 6. BUSCO Analysis of Assembly
#Load BUSCO Module#
#module load GCC/8.3.0 OpenMPI/3.1.4 BUSCO/4.1.2-Python-3.7.4
##Assembly.fasta BUSCO Analysis##
#genome_file='/scratch/user/devonjboland/B-Race/assembly.fasta'
#busco_mode='genome'                 # genome, transcriptome, proteins
#busco_lineage='/scratch/user/devonjboland/odb10/lineages/chlorophyta_odb10'
#lineage_name=$(basename $busco_lineage)
augustus_species='chlamydomonas'
augustus_model='partial'            # partial (default), intronless, complete, atleastone, exactlyone
out_prefix="out_busco4_${lineage_name}_${augustus_species}.Assembly"
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
## Pilonx1.fasta BUSCO Analysis##
genome_file='/scratch/user/devonjboland/B-Race/assembly-pilonx1.fasta'
busco_mode='genome'                 # genome, transcriptome, proteins
busco_lineage='/scratch/user/devonjboland/odb10/lineages/chlorophyta_odb10'
lineage_name=$(basename $busco_lineage)
augustus_species='chlamydomonas'
augustus_model='partial'            # partial (default), intronless, complete, atleastone, exactlyone
out_prefix="out_busco4_${lineage_name}_${augustus_species}.Pilonx1"
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
## Pilonx2.fasta BUSCO Analysis##
genome_file='/scratch/user/devonjboland/B-Race/assembly-pilonx2.fasta'
busco_mode='genome'                 # genome, transcriptome, proteins
busco_lineage='/scratch/user/devonjboland/odb10/lineages/chlorophyta_odb10'
lineage_name=$(basename $busco_lineage)
augustus_species='chlamydomonas'
augustus_model='partial'            # partial (default), intronless, complete, atleastone, exactlyone
out_prefix="out_busco4_${lineage_name}_${augustus_species}.Pilonx2"
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
 