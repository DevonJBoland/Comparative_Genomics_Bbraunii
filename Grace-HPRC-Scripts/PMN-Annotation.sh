#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=PMN-Annotation-Test    #Set the job name to "JobExample2"
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

## Program Variables ##
e2p2='/sw/eb/sw/E2P2/3.1-Java-1.8/runE2P2.v3.1.py'
proteins='Arace_EST_Proteome_v1.fasta'
## Load E2P2 ##
ml E2P2/3.1-Java-1.8
## Run E2P2 Pipeline on Arace Protein Sequences ##
mkdir E2P2
python ${e2p2} -i ${proteins} -o E2P2

## Load SAVI ##
#ml SAVI/3.1-Java-1.8
## Load PerlCyc ##
#ml GCCcore/9.3.0 PerlCyc/1.21-Perl-5.30.2
## Load PlantClusterFinder ##
#ml PlantClusterFinder/1.0-Java-1.8
## Load JavaCyc ##
#ml JavaCyc/1.0-Java-1.
