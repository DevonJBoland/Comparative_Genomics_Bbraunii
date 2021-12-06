# You have installed both python3 and the required python  packages

#  To install the required packaged
# cd to the downloaded folder and run:
pip install -r requiremennts.txt

# To test to validate if the program will work correctly
# cd to Test directory in the program folder and run the shell script with bash
# This shell script will run blast between rna.fa and the genome  (TAIR10.1).

# The script itself will then exectute the python extraction script
# that will pull out the aligned sequence plus upto 1000 bases flanking either side
# if possible.

# The python script will output two files, the fasta file containing the extracted
# sequences, and a log file sumarazing the length of extracted seqeunnces, and the
# original blastn alignment length.