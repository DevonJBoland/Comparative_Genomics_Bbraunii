

# Load Hisat2
ml HISAT2/2.2.1-foss-2018b-Python-3.6.6

# Index Assembly

# Run Alignment

# Unload Hisat2
ml purge

# Load SAMtools
ml SAMtools/1.9-intel-2018b

# Convert Sam to Bam File

# Remove large uneeded sam file
rm
# Sort Bam file

# Index bam file

# Unload samtools
ml purge

# Load BESST
ml BESST/2.0-ictce-7.1.2-Python-2.7.6

