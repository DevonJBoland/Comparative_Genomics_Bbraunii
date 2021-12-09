#!/Users/devonboland/anaconda3/envs/Python-Bubbles/bin/python3

import pandas as pd
import re
import concurrent.futures as cf



# Function to format GTF data by feature
# Need to add a filtering step where genes are removed if they do not contain all features
def screen(x):
    col_names = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    data = pd.read_csv(x, index_col=False, header=None, delimiter="\t", names=col_names)
    # Prefiltering Step
    incomplete_genes_df = pd.DataFrame(columns=col_names)
    complete_genes_df = pd.DataFrame(columns=col_names)
    transcript_ID = [data.at[i, 'attribute'] for i in data.index if data.at[i, 'feature'] == 'transcript']

    def threaded_screen(transcript):
        temp_df = pd.DataFrame(columns=col_names)
        for y in range(0, len(data.index)):
            if re.search(transcript, data.at[y, 'attribute']):
                entry = data.loc[y]
                temp_df = temp_df.append(entry, ignore_index=True)
        feature_checklist = [temp_df.at[i, 'feature'] for i in temp_df.index]
        if feature_checklist.count('5\'-UTR') >= 1 and feature_checklist.count('start_codon') == 1 and \
                feature_checklist.count('CDS') >= 1 and feature_checklist.count('stop_codon') == 1 and \
                feature_checklist.count('3\'-UTR') >= 1:
            complete_genes_df = complete_genes_df.append(temp_df, ignore_index=True)
        else:
            incomplete_genes_df = incomplete_genes_df.append(temp_df, ignore_index=True)


    with cf.ProcessPoolExecutor() as executor:
        executor.map(threaded_screen, transcript_ID)

    complete_genes_df.to_csv("Arace-complete-filtered.gtf", sep="\t", header=True, index=False)
    incomplete_genes_df.to_csv("Arace-incomplete-filtered.gtf", sep="\t", header=True, index=False)

 """   # Extraction Step
    
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
        elif data.at[i, 2] == '5\'-UTR':

        elif data.at[i, 2] == 'start_codon':

        elif data.at[i, 2] == 'CDS':

        elif data.at[i, 2] == 'stop_codon':

        elif data.at[i, 2] == '3\'-UTR':
"""



with open("UTR-GTFs/Arace-augustus.hints_utr.gtf", "r") as file:
    screen(file)
