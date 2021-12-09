#!/Users/devonboland/anaconda3/envs/Python-Bubbles/bin/python3

import pandas as pd
import re
import concurrent.futures as cf

# Function to format GTF data by feature
# Need to add a filtering step where genes are removed if they do not contain all features


def format_data(x):
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

        if isnegeativegene(temp_df):
            temp_df = reversechunck(temp_df)
            print(temp_df)




def isnegeativegene(gene_chunck):
    return gene_chunck.at[0, 'strand'] == '-'


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
    format_data(file)
