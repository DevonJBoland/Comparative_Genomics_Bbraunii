#!/bin/bash

# Wrapper script to create intron length and frequency from gtf files

# Assign Variables
input_directory=($1)
output_directory=($2)
# These you may need to change
venv='/Users/devonboland/Hybrid-Assembly/PythonScripts/venv'
Extract_IntronGFF_Count='/Users/devonboland/Hybrid-Assembly/PythonScripts/Intron_Length_Analysis/Extract_IntronGFF_Count.py'
# Activate python environment to execute extraction script
source ${venv}/bin/activate
# Create list of files in input directory to loop through later
mkdir $2
ls $1 > loop-through.txt
cat loop-through.txt | while read line
do
python3 ${Extract_IntronGFF_Count} ${1}${line} ${2}${line}-count.csv
done
rm loop-through.txt
deactivate
