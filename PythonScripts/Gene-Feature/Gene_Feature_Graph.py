#!/Users/devonboland/anaconda3/envs/Python-Bubbles/bin/python3

import pandas as pd
import matplotlib.pyplot as plt
import re
import argparse as ap

col_names = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

def import_data(data):
    data_array = pd.read_csv(data, delimiter='\t', index_col=False)
    return data_array


def chunck_gene(annotation_df):
    transcript_ID = [annotation_df.at[n, 'attribute'] for n in annotation_df.index
                     if annotation_df.at[n, 'feature'] == 'transcript']

    df_dict = {}
    start = 0
    # for transcript in range(len(transcript_ID) - 1):
    for transcript in range(0 ,2, 1): # for testing range max = 2
        df = pd.DataFrame(columns=col_names)
        # for n in range(start, len(annotation_df)-1):
        for n in range(start, 72): # for testing range max = 72
            if re.search(transcript_ID[transcript], annotation_df.at[n, 'attribute']):
                entry = annotation_df.loc[n]
                df = df.append(entry, ignore_index=True)
            else:
                start = n
                break
        df_dict[transcript_ID[transcript]] = df

    return df_dict, transcript_ID


def intron_distribution(variable): # rename the input vairable to be relevant




with open("UTR-GTFs/Lrace-augustus.hints_utr.gtf-intron", 'r') as file:
    annotations = import_data(file)
    transcripts, list_transcript = chunck_gene(annotations)
