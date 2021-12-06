#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=InterProScan-GO-Annotation    #Set the job name to "JobExample2"
#SBATCH --time=01-00:00:00           # max job run time dd-hh:mm:ss
#SBATCH --nodes=1
##SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=40          # CPUs (threads) per command
#SBATCH --mem=2900G                  # total memory per node
#SBATCH --partition=bigmem
#SBATCH --output=stdout.%j.txt          # save stdout to file
#SBATCH --error=stderr.%j.txt           # save stderr to file

##SET OPTIONAL DIRECTIVES
#SBATCH --mail-type=ALL              #Send email on all job events
#SBATCH --mail-user=devonjboland@tamu.edu    #Send all emails to

# Assign Variables
Arace='/scratch/user/devonjboland/GO-Annotation/BobraA-cor.fasta'
Brace='/scratch/user/devonjboland/GO-Annotation/Bbraunii_502_v2.1.protein-cor.fa'
Lrace='/scratch/user/devonjboland/GO-Annotation/BobraL-cor.fasta'
#threads=$SLURM_CPUS_PER_TASK
#interpro=$EBROOTINTERPROSCAN
#temp='/scratch/user/devonjboland/GO-Annotation/'
applications='TIGRFAM,PIRSF,SMART,PrositeProfiles,PrositePatterns,HAMAP,PfamA,PRINTS,SuperFamily,Coils,SignalP-EUK,TMHMM'

# Load InterProScan
ml GCC/9.3.0  OpenMPI/4.0.3 InterProScan/5.52-86.0
#ml GCC/9.3.0  OpenMPI/4.0.3  InterProScan/5.51-85.0

# Command to execute interproscan
# java -Xmx2900G -jar ${interpo}interproscan-5.jar -cpu ${threads} -f XML  -goterms -i  -iprlookup -o  -pa -T ${temp}

#tamulauncher --commands-pernode 8 interpro-commands.txt

# Arace
interproscan.sh --disable-precalc --tempdir $TMPDIR --input ${Arace} --applications $applications -o Arace --formats XML -cpu 40
# Brace
interproscan.sh --disable-precalc --tempdir $TMPDIR --input ${Brace} --applications $applications -o Brace --formats XML -cpu 40
# Lrace
interproscan.sh --disable-precalc --tempdir $TMPDIR --input ${Lrace} --applications $applications -o Lrace --formats XML -cpu 40