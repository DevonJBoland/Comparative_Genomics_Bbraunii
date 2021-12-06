#!/bin/bash
<<README
    - BUSCO manual: https://busco.ezlab.org/busco_userguide.html
    - Homepage: http://busco.ezlab.org
README

################################### VARIABLES ##################################
# TODO Edit these variables as needed:

########## INPUTS ##########
genome_file='/scratch/user/devonjboland/Canu-3TBNode-CanuAssembly/Pilon-ContigPolishing/Merged-pilon/ Merged+Canu_PILONx2.scaffolds.fasta'

######## PARAMETERS ########
threads=80                          # make sure this is <= your BSUB -n value

# see available BUSCO 4.1.2 lineages in this directory: /scratch/data/bio/BUSCO/odb10/lineages/
busco_mode='genome'                 # genome, transcriptome, proteins
busco_lineage='/scratch/user/devonjboland/Canu-3TBNode-CanuAssembly/BUSCO/download_lineage/lineages/chlorophyta_odb10'
lineage_name=$(basename $busco_lineage)

# --species:  see available species here: /sw/eb/sw/AUGUSTUS/3.3.3-foss-2019b/config/species/
augustus_species='chlamydomonas'
augustus_model='partial'            # partial (default), intronless, complete, atleastone, exactlyone

########## OUTPUTS #########
out_prefix="out_busco4_${lineage_name}_${augustus_species}.Merged-PILONx2"

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

busco --offline --in $genome_file --mode $busco_mode --augustus_species $augustus_species --cpu $threads --lineage $busco_lineage --out $out_prefix --augustus_parameters="--genemodel=$augustus_model"

################################################################################
<<CITATION
    - Acknowledge TAMU HPRC: https://hprc.tamu.edu/research/citations.html

    - BUSCO:
        BUSCO: assessing genome assembly and annotation completeness with single-copy orthologs.
        Felipe A. Simão, Robert M. Waterhouse, Panagiotis Ioannidis, Evgenia V. Kriventseva, and Evgeny M. Zdobnov
        Bioinformatics, published online June 9, 2015 | doi: 10.1093/bioinformatics/btv351
CITATION
