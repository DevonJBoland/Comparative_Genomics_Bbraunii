#!/bin/bash

# Variables to Set Prior to Use
genome='GCF_000001735.4_TAIR10.1_genomic.fna'
query='rna.fna'
extract='/Users/devonboland/Hybrid-Assembly/PythonScripts/TE_Model_Identification.py'
hits='formatted_hits.csv'
output='Final_Reslts.fa'

# Format Genome as Blast Nucleotide Database
makeblastdb -in ${genome} -dbtype nucl

# Blast Query Sequences Against Genome
blastn -query ${query} -db ${genome} -outfmt 6 -out blastn-results.csv
cat blastn-results.csv | cut -f2,9,10 | sed 's/\t/,/g' > ${hits}
# Code up to the previous line is tested and working {Devon J Boland}
# Run Custom Extraction Script to extract target aligned sequences + 1000bp 
# flanking sequence
python3 ${extract} ${genome} ${hits} ${output}
