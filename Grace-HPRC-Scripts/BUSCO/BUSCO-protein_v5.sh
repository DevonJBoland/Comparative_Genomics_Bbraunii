#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=Pytozome-BUSO            # job name
#SBATCH --time=2-00:00:00           # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=80          # CPUs (threads) per command
#SBATCH --mem=2900G                  # total memory per node
#SBATCH --partition=bigmem
#SBATCH --output=stdout.%j.txt          # save stdout to file
#SBATCH --error=stderr.%j.txt           # save stderr to file

module load GCC/10.2.0 OpenMPI/4.0.5 BUSCO/5.1.2

cat names.txt | while read line
do
################################### VARIABLES ##################################
# TODO Edit these variables as needed:
########## INPUTS ##########
genome_file=${line}
######## PARAMETERS ########
threads=$SLURM_CPUS_PER_TASK
busco_mode='proteins'                 # genome, transcriptome, proteins
busco_lineage='/scratch/user/devonjboland/odb10/lineages/chlorophyta_odb10'         # check available databases at /scratch/data/bio/busco5/lineages/
#download_path='/scratch/data/bio/busco5'    # use this if database is available
#augustus_species='chlamydomonas' # see available species here: /sw/eb/sw/AUGUSTUS/3.4.0-foss-2020b/config/species/
#augustus_model='partial'            # partial (default), intronless, complete, atleastone, exactlyone
########## OUTPUTS #########
out_prefix=${line}-output
########## CONFIG ##########
export AUGUSTUS_CONFIG_PATH="$SCRATCH/my_augustus_config/config"
################################### COMMANDS ###################################
busco --offline --in $genome_file --mode $busco_mode --cpu $threads --lineage $busco_lineage --out $out_prefix
################################################################################
wait
done
