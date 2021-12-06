"""
Purpose: To extract and determine the count(frequency) of intron length present from gtf files.
Created: 08-11-21
Created by: Devon J. Boland
"""
import pandas as pd
import csv as csv
import argparse as ap

# Format input from user using argparser
parser = ap.ArgumentParser(description="Program to extract intron length and frequencies")
parser.add_argument('input', type=str, help="Genome feature file in gtf format (ex:gene.gtf)")
parser.add_argument('feature', type=str, help="Genome feature to extract (ex: gene, intron, exon, ")
parser.add_argument("output", type=str, help="Output file name (ex:output.csv)")
args = parser.parse_args()

# Open the gtf file and import gtf contents into a pandas dataframe
with open(args.input, "r") as gtf:
    data = pd.read_csv(gtf, sep='\t', index_col=False, header=None)
    data_tuple = [tuple(row) for row in data.values] # enumerate through the df and create a tuple with a dictionary for each row

# Extract all lines from tuple that correspond to an intron and keep the 0, 2, 3, 4 columns
limit = len(data_tuple)
name = []
state = []
length = []
for i in range(0, limit):
    if str(data_tuple[i][7]) == str('intron'):
        name.append(str(data_tuple[i][0]))
        math = int(data_tuple[i][4]) - int(data_tuple[i][3])
        length.append(str(math))

# Create a pandas data frame from the collected intron data
introns_df = pd.DataFrame(list(zip(name, length)), columns=['Unitig Location', 'Feature Length'])
# Count frequency of intron length and cast to frequency dict then sorty by key converting to a list and write to csv file
with open(args.output, "w") as log:
    frequency = introns_df['Feature Length'].value_counts().to_dict()
    int_frequency_items = {int(k):int(v) for k, v in frequency.items()}
    sorted_frequency = sorted(int_frequency_items.items())
    writer = csv.writer(log)
    for row in sorted_frequency:
        writer.writerow(row)
# Currently works on test data