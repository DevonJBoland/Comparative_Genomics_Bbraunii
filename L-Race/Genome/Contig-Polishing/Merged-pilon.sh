#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=PilonPolishx2      #Set the job name to "JobExample2"
#SBATCH --time=2-0:00:00               #Set the wall clock limit to 6hr and 30min
#SBATCH --ntasks=80
#SBATCH --ntasks-per-node=80          #Request 20 tasks/cores per node
#SBATCH --mem=2900G                     #Request 2GB per node
#SBATCH --partition=bigmem
#SBATCH --output=stdout.%j.txt          # save stdout to file
#SBATCH --error=stderr.%j.txt           # save stderr to file 

##SET OPTIONAL DIRECTIVES
#SBATCH --mail-type=END              #Send email on all job events
#SBATCH --mail-user=devonjboland@tamu.edu    #Send all emails to email_address

# Set Variables and Directories
# Direcotry variables
workdir=/scratch/user/devonjboland/Canu-3TBNode-CanuAssembly
outputdir=$workdir/Pilon-ContigPolishing/Merged-pilon
# Illumina paths
illuminareads=$workdir/IlluminaDataSet
pe1=$illuminareads/Race_L_R1.fastq.gz
pe2=$illuminareads/Race_L_R2.fastq.gz

# Load Bowtie2 Module
ml GCC/7.3.0-2.30
ml OpenMPI/3.1.1
ml HISAT2/2.2.0
ml SAMtools/1.9
ml Pilon/1.23-Java-11

########
# First round of mapping and polishing
########
# First round of polishing variables
contigs=$workdir/Pilon-ContigPolishing/Merged-pilon
ht2index=Mergedcontigs
frags=$workdir/Pilon-ContigPolishing/Merged-pilon/IlluminaAligned-Sorted.bam
output=Merged+Canu_PILON.scaffolds

# Build HISAT2 Index Files of L race contigs
hisat2-build -p 80 $contigs $ht2index
wait

# Map reads to contigs
hisat2  --sensitive -p 80 -x $ht2index -1 $pe1 -2 $pe2 -S IlluminaAligned.sam
samtools view -@ 80 -bSh -o IlluminaAligned.bam IlluminaAligned.sam

# Remove Large Sam file
#rm IlluminaAligned.sam

# Sort BAM file
samtools sort -o IlluminaAligned-Sorted -O bam -@ 80 IlluminaAligned.bam
mv IlluminaAligned-Sorted IlluminaAligned-Sorted.bam

# Index Bam File
samtools index -b -@ 80 $frags $frags.bai

# Correct ONT Contigs with Pilon Using Illumina Mapped Reads
java -Xmx2900G -jar $EBROOTPILON/pilon.jar --threads 80 --genome $contigs --frags $frags --output $output --outdir $outputdir

########
# Second round of mapping and polishing
########
# Second round of polishing variables
contigsPilon=$workdir/Pilon-ContigPolishing/Merged-pilon/Merged+Canu_PILON.scaffolds.fasta
ht2index2=Mergedcontigs2
frags2=$workdir/Pilon-ContigPolishing/Merged-pilon/IlluminaAligned-Sorted2.bam
output2=Merged+Canu_PILONx2.scaffolds

# Build HISAT2 Index Files of L race contigs
hisat2-build -p 80 $contigsPilon $ht2index2
wait

# Map reads to contigs
hisat2  --sensitive -p 80 -x $ht2index2 -1 $pe1 -2 $pe2 -S IlluminaAligned2.sam
samtools view -@ 80 -bSh -o IlluminaAligned2.bam IlluminaAligned2.sam

# Remove Large Sam file
#rm IlluminaAligned2.sam

# Sort BAM file
samtools sort -o IlluminaAligned-Sorted2 -O bam -@ 80 IlluminaAligned2.bam
mv IlluminaAligned-Sorted2 IlluminaAligned-Sorted2.bam

# Index Bam File
samtools index -b -@ 80 $frags2 $frags2.bai

# Correct ONT Contigs with Pilon Using Illumina Mapped Reads
java -Xmx2900G -jar $EBROOTPILON/pilon.jar --threads 80 --genome $contigsPilon --frags $frags2 --output $output2 --outdir $outputdir
