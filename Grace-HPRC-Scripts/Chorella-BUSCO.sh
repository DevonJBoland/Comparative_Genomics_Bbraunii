#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=Chlorella-BUSCO            # job name
#SBATCH --time=2-00:00:00           # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=80          # CPUs (threads) per command
#SBATCH --mem=2900G                  # total memory per node
#SBATCH --partition=bigmem
#SBATCH --output=stdout.%j.txt          # save stdout to file
#SBATCH --error=stderr.%j.txt           # save stderr to file

##SET OPTIONAL DIRECTIVES
#SBATCH --mail-type=ALL              #Send email on all job events
#SBATCH --mail-user=devonjboland@tamu.edu    #Send all emails to

module load GCC/10.2.0 OpenMPI/4.0.5 BUSCO/5.1.2
################################### VARIABLES ##################################
# TODO Edit these variables as needed:
########## INPUTS ##########
genome_file='ChloA99_1_AssemblyScaffolds_Repeatmasked.fasta'
######## PARAMETERS ########
threads=$SLURM_CPUS_PER_TASK
busco_mode='genome'                 # genome, transcriptome, proteins
busco_lineage='/scratch/user/devonjboland/odb10/lineages/chlorophyta_odb10'         # check available databases at /scratch/data/bio/busco5/lineages/
#download_path='/scratch/data/bio/busco5'    # use this if database is available
augustus_species='chlamydomonas' # see available species here: /sw/eb/sw/AUGUSTUS/3.4.0-foss-2020b/config/species/
augustus_model='partial'            # partial (default), intronless, complete, atleastone, exactlyone
########## OUTPUTS #########
out_prefix="out_busco5_chlorophyta_odb10_chlamydomonas.Chlorella"
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
