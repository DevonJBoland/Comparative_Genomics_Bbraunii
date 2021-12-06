#!/Users/devonboland/anaconda3/envs/Python-Bubbles/bin/python3

import pandas as pd

# Function to format GTF data by feature
def extract(x):
    data = pd.read_csv(x, index_col=False, header=None, delimiter="\t")
    #print(data[2])
    col_names = ['Unitig','Feature','strand','']
    gene = pd.DataFrame()
    five_UTR = pd.DataFrame()
    start = pd.DataFrame()
    CDS = pd.DataFrame()
    stop = pd.DataFrame()
    three_UTR = pd.DataFrame()
    s
    for i in range(0,len(data.index)):
        if data.at[i, 2] == 'gene':
            new_row = data.loc[i]
            gene.append(new_row, ignore_index=True)
        elif data.at[i, 2] == '5\'-UTR':

        elif data.at[i, 2] == 'start_codon':

        elif data.at[i, 2] == 'CDS':

        elif data.at[i, 2] == 'stop_codon':

        elif data.at[i, 2] == '3\'-UTR':


    print(gene)

with open("UTR-GTFs/Arace-augustus.hints_utr.gtf", "r") as file:
    extract(file)




