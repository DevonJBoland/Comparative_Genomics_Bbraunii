#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=BRAKER-UTR    #Set the job name to "JobExample2"
#SBATCH --time=02-00:00:00           # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=48          # CPUs (threads) per command
#SBATCH --mem=2900G                  # total memory per node
#SBATCH --partition=bigmem
#SBATCH --output=stdout.%j.txt          # save stdout to file
#SBATCH --error=stderr.%j.txt           # save stderr to file

##SET OPTIONAL DIRECTIVES
#SBATCH --mail-type=ALL              #Send email on all job events
#SBATCH --mail-user=devonjboland@tamu.edu    #Send all emails to

# Assign Variables
chlamy_genome="Chlre5_6_AssemblyScaffolds_Repeatmasked.fasta"
chlamy_bam="mapped-trim-reads.bam"
chlamy_workdir="/scratch/user/devonjboland/Annotation/Creinhardtii"

chlorella_genome="ChlNC64A_1_nuclear_scaffolds_softmasked.fasta"
chlorella_bam="mapped-trim-reads.bam"
chlorella_workdir="/scratch/user/devonjboland/Annotation/Cvariabilis"

coccomyxa_genome="Coccomyxa_C169_v2_masked_genomic_scaffolds.fasta"
coccomyxa_bam="mapped-trim-reads.bam"
coccomyxa_workdir="/scratch/user/devonjboland/Annotation/Csubellipsoidea"

augustus_path="$SCRATCH/my_augustus_config/config"
threads=$SLURM_CPUS_PER_TASK

# Load BRAKER Modules
ml GCC/9.3.0  OpenMPI/4.0.3  BRAKER/2.1.6-Python-3.8.2
export AUGUSTUS_CONFIG_PATH="$SCRATCH/my_augustus_config/config"
export GUSHR_PATH=$SCRATCH/GUSHR

# Run BRAKER with Mapped RNA-Seq Data Creinhardtii
#cd ${chlamy_workdir}

#braker.pl --GUSHR_PATH=$SCRATCH/GUSHR --AUGUSTUS_CONFIG_PATH=${augustus_path} --softmasking --addUTR=on --species=Arace --genome=${Aracegenome} --bam=${Aracernaseq_bam} --cores=${threads} # Commented out as it did not work as well with similar code blocks.

#braker.pl --GUSHR_PATH=$SCRATCH/GUSHR --AUGUSTUS_CONFIG_PATH=${augustus_path} --genome=${chlamy_genome} --addUTR=on --softmasking \
    --bam=${chlamy_bam} --workingdir=${chlamy_workdir} --cores=8 \
    --species=Creinhardtii

# Run BRAKER with Mapped RNA-Seq Data Cvariabilis
cd ${chlorella_workdir}
#braker.pl --GUSHR_PATH=$SCRATCH/GUSHR --AUGUSTUS_CONFIG_PATH=${augustus_path} --softmasking --addUTR=on --species=Lrace --genome=${Lracegenome} --bam=${Lracernaseq_bam} --cores=${threads}

braker.pl --GUSHR_PATH=$SCRATCH/GUSHR --AUGUSTUS_CONFIG_PATH=${augustus_path} --genome=${chlorella_genome} --addUTR=on --softmasking \
    --bam=${chlorella_bam} --workingdir=${chlorella_workdir} --cores=8 \
    --species=Cvariabilis

# Run BRAKER with Mapped RNA-Seq Data Csubellipsoidea
cd ${coccomyxa_workdir}
#braker.pl --GUSHR_PATH=$SCRATCH/GUSHR --AUGUSTUS_CONFIG_PATH=${augustus_path} --softmasking --addUTR=on --species=Brace --genome=${Bracegenome} --bam=${Bracernaseq_bam} --cores=${threads}

braker.pl --GUSHR_PATH=$SCRATCH/GUSHR --AUGUSTUS_CONFIG_PATH=${augustus_path} --genome=${coccomyxa_genome} --addUTR=on --softmasking \
    --bam=${coccomyxa_bam} --workingdir=${coccomyxa_workdir} --cores=8 \
    --species=Csubellipsoidea
