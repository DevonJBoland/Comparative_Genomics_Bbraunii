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
chlamy_genome=""
chlamy_bam=""
chlamy_workdir=""

chlorella_genome=""
chlorella_bam=""
chlorella_workdir=""

coccomyxa_genome=""
coccomyxa_bam=""
coccomyxa_workdir=""

augustus_path="$SCRATCH/my_augustus_config/config"
threads=$SLURM_CPUS_PER_TASK
# Load BRAKER Modules
ml GCC/9.3.0  OpenMPI/4.0.3  BRAKER/2.1.6-Python-3.8.2
export AUGUSTUS_CONFIG_PATH="$SCRATCH/my_augustus_config/config"
export GUSHR_PATH=$SCRATCH/GUSHR

# Run BRAKER with Mapped RNA-Seq Data Arace
cd ${Araceworkdir}

#braker.pl --GUSHR_PATH=$SCRATCH/GUSHR --AUGUSTUS_CONFIG_PATH=${augustus_path} --softmasking --addUTR=on --species=Arace --genome=${Aracegenome} --bam=${Aracernaseq_bam} --cores=${threads} # Commented out as it did not work as well with similar code blocks.

braker.pl --GUSHR_PATH=$SCRATCH/GUSHR --AUGUSTUS_CONFIG_PATH=${augustus_path} --genome=${Aracegenome} --addUTR=on --softmasking \
    --bam=${Aracernaseq_bam} --workingdir=${Araceworkdir} \
    --AUGUSTUS_hints_preds=${Araceworkdir}/braker/augustus.hints.gtf --cores=8 \
    --skipAllTraining --species=Arace

# Run BRAKER with Mapped RNA-Seq Data Lrace
cd ${Lraceworkdir}
#braker.pl --GUSHR_PATH=$SCRATCH/GUSHR --AUGUSTUS_CONFIG_PATH=${augustus_path} --softmasking --addUTR=on --species=Lrace --genome=${Lracegenome} --bam=${Lracernaseq_bam} --cores=${threads}

braker.pl --GUSHR_PATH=$SCRATCH/GUSHR --AUGUSTUS_CONFIG_PATH=${augustus_path} --genome=${Lracegenome} --addUTR=on --softmasking \
    --bam=${Lracernaseq_bam} --workingdir=${Lraceworkdir} \
    --AUGUSTUS_hints_preds=${Lraceworkdir}/braker/augustus.hints.gtf --cores=8 \
    --skipAllTraining --species=Lrace

# Run BRAKER with Mapped RNA-Seq Data Lrace
cd ${Braceworkdir}
#braker.pl --GUSHR_PATH=$SCRATCH/GUSHR --AUGUSTUS_CONFIG_PATH=${augustus_path} --softmasking --addUTR=on --species=Brace --genome=${Bracegenome} --bam=${Bracernaseq_bam} --cores=${threads}

braker.pl --GUSHR_PATH=$SCRATCH/GUSHR --AUGUSTUS_CONFIG_PATH=${augustus_path} --genome=${Bracegenome} --addUTR=on --softmasking \
    --bam=${Bracernaseq_bam} --workingdir=${Braceworkdir} \
    --AUGUSTUS_hints_preds=${Braceworkdir}/braker/augustus.hints.gtf --cores=8 \
    --skipAllTraining --species=Brace
