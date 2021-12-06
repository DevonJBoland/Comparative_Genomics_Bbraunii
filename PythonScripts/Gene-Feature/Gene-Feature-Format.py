#!/Users/devonboland/anaconda3/envs/Python-Bubbles/bin/python3

import pandas as pd

# Function to format GTF data by feature
# Need to add a filtering step where genes are removed if they do not contain all features
def extract(x):
    col_names = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    data = pd.read_csv(x, index_col=False, header=None, delimiter="\t", names=col_names)
    gene = pd.DataFrame(columns=col_names)
    five_UTR = pd.DataFrame(columns=col_names)
    start = pd.DataFrame(columns=col_names)
    CDS = pd.DataFrame(columns=col_names)
    stop = pd.DataFrame(columns=col_names)
    three_UTR = pd.DataFrame(columns=col_names)

    for i in range(0,len(data.index)): # Still need to add decision on which order to add or substract
        if data.at[i, 'feature'] == 'gene':
            entry = data.loc[i]
            gene = gene.append(entry, ignore_index=True)
        """elif data.at[i, 2] == '5\'-UTR':

        elif data.at[i, 2] == 'start_codon':

        elif data.at[i, 2] == 'CDS':

        elif data.at[i, 2] == 'stop_codon':

        elif data.at[i, 2] == '3\'-UTR':
"""


with open("UTR-GTFs/Arace-augustus.hints_utr.gtf", "r") as file:
    extract(file)
