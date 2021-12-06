#!/bin/bash

# Small script to combine tree files from IQTREE into a single file.
# Remove low supported branches, and create a single alternative species tree

# Environment was created in conda

# Combine tree files into a single tree file
cat *.treefile > chlorophyta_algae_1189.gene.tree

# Use newick_utilites to remove low supported branches
nw_ed  chlorophyta_algae_1189.gene.tree 'i & b<=10' o > chlorophyta_algae_BS10.gene.tree

# Run ASTRAL-III to combine into a single tree
java -jar astral.5.7.7.jar -i chlorophyta_algae_BS10.gene.tre -o chlorophyta_BS10.gene.tre 2> chlorophyta_bs10.gene.tre.log
