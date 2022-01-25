#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=alphafold2_CYP153A        # job name
#SBATCH --time=1-00:00:00           # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=24          # CPUs (threads) per command
#SBATCH --mem=180G                  # total memory per node
#SBATCH --gres=gpu:a100:1           # request 1 A100 GPU
#SBATCH --output=stdout.%j.txt          # save stdout to file
#SBATCH --error=stderr.%j.txt           # save stderr to file

##SET OPTIONAL DIRECTIVES
#SBATCH --mail-type=ALL              #Send email on all job events
#SBATCH --mail-user=devonjboland@tamu.edu    #Send all emails to

# Assign Variables
DOWNLOAD_DIR=/scratch/data/bio/alphafold
fasta="CYP153A6.fasta"
outdir="/scratch/user/devonjboland/Dakota_Grote_CYP153A6"


# Run Alphafold2 on Input Fasta
singularity exec --nv /sw/hprc/sw/bio/containers/alphafold_2.1.1.sif python /app/alphafold/run_alphafold.py  \
  --fasta_paths=${fasta}  \
  --data_dir=$DOWNLOAD_DIR  \
  --uniref90_database_path=$DOWNLOAD_DIR/uniref90/uniref90.fasta  \
  --mgnify_database_path=$DOWNLOAD_DIR/mgnfy/mgy_clusters_2018_12.fa  \
  --bfd_database_path=$DOWNLOAD_DIR/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt  \
  --uniclust30_database_path=$DOWNLOAD_DIR/uniclust30/uniclust30_2018_08/uniclust30_2018_08  \
  --pdb70_database_path=$DOWNLOAD_DIR/pdb70/pdb70  \
  --template_mmcif_dir=$DOWNLOAD_DIR/pdb_mmcif/mmcif_files  \
  --obsolete_pdbs_path=$DOWNLOAD_DIR/pdb_mmcif/obsolete.dat  \
  --output_dir=${outdir}  \
  --model_preset=monomer  \
  --db_preset=full_dbs  \
  --max_template_date=2022-01-01