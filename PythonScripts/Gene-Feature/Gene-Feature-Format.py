#!/Users/devonboland/anaconda3/envs/Python-Bubbles/bin/python3

import pandas as pd
import re
from toolz import interleave
#import concurrent.futures as cf

# Function to format GTF data by feature
# Need to add a filtering step where genes are removed if they do not contain all features


def filter_genes(x):
    col_names = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    data = pd.read_csv(x, index_col=False, header=None, delimiter="\t", names=col_names)
    # Prefiltering Step
    incomplete_genes_df = pd.DataFrame(columns=col_names)
    complete_genes_df = pd.DataFrame(columns=col_names)
    transcript_ID = [data.at[i, 'attribute'] for i in data.index if data.at[i, 'feature'] == 'transcript']

    for transcript in transcript_ID:
        temp_df = pd.DataFrame(columns=col_names)
        for y in range(0, len(data.index)):
            if re.search(transcript, data.at[y, 'attribute']):
                entry = data.loc[y]
                temp_df = temp_df.append(entry, ignore_index=True)

        feature_checklist = [temp_df.at[i, 'feature'] for i in temp_df.index]
        if meetsgenecriteria(feature_checklist):
            complete_genes_df = complete_genes_df.append(temp_df, ignore_index=True)
        else:
            incomplete_genes_df = incomplete_genes_df.append(temp_df, ignore_index=True)

        to_remove = temp_df.at[1, 'attribute']
        data = data[data.attribute != to_remove]
        if data.at[0, 'feature'] == 'gene':
            data = data.iloc[2:, :]
        else:
            data = data.iloc[1:, :]
        data = data.reset_index(drop=True)

    # complete_genes_df.to_csv("test-complete-filtered.gtf", sep="\t", header=True, index=False)
    # incomplete_genes_df.to_csv("test-incomplete-filtered.gtf", sep="\t", header=True, index=False)

    return complete_genes_df


def get_intron(gene_df):  # Not working properly to create a list of transcript IDs
    col_names = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    transcript_ID = [gene_df.at[n, 'attribute'] for n in gene_df.index if gene_df.at[n, 'feature'] == 'transcript']
    # if not transcript_ID:
    #     print("transcript list is empty")
    #     return
    for transcript in transcript_ID:
        temp_gene_df = pd.DataFrame(columns=col_names)
        for y in range(0, len(gene_df.index)):
            if re.search(transcript, gene_df.at[y, 'attribute']):
                entry = gene_df.loc[y]
                temp_gene_df = temp_gene_df.append(entry, ignore_index=True)
        if isnegeativegene(temp_gene_df):
            temp_df = reversechunck(temp_gene_df)
        else:
            temp_df = temp_gene_df
        fiveprime_array = temp_gene_df[temp_gene_df.feature == '5\'-UTR']
        cds_array = temp_gene_df[temp_gene_df.feature == 'CDS']
        threeprime_array = temp_gene_df[temp_gene_df.feature == '3\'-UTR']
        seqname = temp_df.at[1, 'seqname']
        source = 'Custom_Intron_Script'
        feature = 'intron'
        score = '.'
        strand = temp_df.at[1, 'strand']
        frame = '.'
        attribute = temp_df.at[1, 'strand']
        temp_intron_df = pd.DataFrame(columns=col_names)
        if len(fiveprime_array) > 1:
            for n in fiveprime_array.index:
                intron_start = int(fiveprime_array.at[n, 'end']+1)
                intron_end = int(fiveprime_array.at[n+1, 'start']-1)
                intron_annotations = [seqname, source, feature, intron_start, intron_end, score, strand, frame,
                                      attribute]
                entry = dict(zip(col_names, intron_annotations))
                temp_intron_df = temp_intron_df.append(entry, ignore_index=True)





        # if isnegeativegene(temp_df):
        #     reverse_df = reversechunck(temp_df)
        #     get_cds_intron(reverse_df)
        # # test = temp_df
        # if isnegeativegene(temp_df):
        #     reverse_temp_df = reversechunck(temp_df)
        #     if meetsgenecriteria(feature_checklist):
        #         return 'Is Complete'
        #     else:
        #         return 'Not Complete'

    

def isnegeativegene(gene_chunck):
    return gene_chunck.at [1, 'strand'] in ['-', '- '] # there is irradic spacing in the strand column throughout the gtf files.

def reversechunck(gene_chunck):
    transcript_array = gene_chunck[gene_chunck.feature == 'transcript']
    features_array = gene_chunck[gene_chunck.feature != 'transcript']
    features_array = features_array.iloc[::-1]
    return transcript_array.append(features_array, ignore_index=True)
    
    


        # if feature_checklist.count('5\'-UTR') >= 1 and feature_checklist.count('start_codon') == 1 and \
        #         feature_checklist.count('CDS') >= 1 and feature_checklist.count('stop_codon') == 1 and \
        #         feature_checklist.count('3\'-UTR') >= 1:
        #     complete_genes_df = complete_genes_df.append(temp_df, ignore_index=True)
        # else:
        #     incomplete_genes_df = incomplete_genes_df.append(temp_df, ignore_index=True)


def meetsgenecriteria(featurelist):
    return featurelist.count('5\'-UTR') >= 1 and featurelist.count('start_codon') == 1 and \
            featurelist.count('CDS') >= 1 and featurelist.count('stop_codon') == 1 and \
            featurelist.count('3\'-UTR') >= 1

    #with cf.ProcessPoolExecutor() as executor:
    #    executor.map(threaded_screen, transcript_ID)

    # complete_genes_df.to_csv("test-complete-filtered.gtf", sep="\t", header=True, index=False)
    # incomplete_genes_df.to_csv("test-incomplete-filtered.gtf", sep="\t", header=True, index=False)


with open("UTR-GTFs/test.gtf", "r") as file:
    completegenes = filter_genes(file)
    # Test if the filtering step was a success
    if completegenes.empty:
        print("Genes Filtered = Failure!")
    else:
        print("Genes Filtered = Success!"+"\n"+f"Output is saved to...")  # Need to insert a string for saved file

    get_intron(completegenes)

