#!/Users/devonboland/anaconda3/envs/Python-Bubbles/bin/python3

import pandas as pd
import re
import argparse as ap

# parser = ap.ArgumentParser(description="Script to filter GTF files off incomplete gene annotations, "
#                                        "and annotate intronic regions")
# parser.add_argument("input", type=str, help="gtf file to extract feature from")
# parser.add_argument("filtered_out", type=str, help="Path to write filtered complete models to")
# parser.add_argument("incomplete_out", type=str, help="Path to write filtered incomplete models to")
# parser.add_argument("intron_annotated_out", type=str, help="Path to write filtered intron annotated "
#                                                               "complete models to")
# args = parser.parse_args()

col_names = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']


def filter_genes(x):
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

    # complete_genes_df.to_csv(args.filtered_out, sep="\t", header=True, index=False)
    # incomplete_genes_df.to_csv(args.incomplete_out, sep="\t", header=True, index=False)

    return complete_genes_df


def get_intron(gene_df):  # Not working properly to create a list of transcript IDs
    transcript_ID = [gene_df.at[n, 'attribute'] for n in gene_df.index if gene_df.at[n, 'feature'] == 'transcript']
    transcript_with_introns_temp = pd.DataFrame(columns=col_names)
    transcript_with_introns = pd.DataFrame(columns=col_names)
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

        # Assign individual arrays for each feature, reset index for concat later
        """This needs to be fixed to handle cases where a UTR may not be present,
        as of now it causes an error when this occurs"""
        transcript_array = temp_df[temp_df.feature == 'transcript']
        transcript_array.reset_index(inplace=True, drop=True)
        fiveprime_array = temp_df[temp_df.feature == '5\'-UTR']
        fiveprime_array.reset_index(inplace=True, drop=True)
        startcodon_array = temp_df[temp_df.feature == 'start_codon']
        startcodon_array.reset_index(inplace=True, drop=True)
        cds_array = temp_df[temp_df.feature == 'CDS']
        cds_array.reset_index(inplace=True, drop=True)
        stopcodon_array = temp_df[temp_df.feature == 'stop_codon']
        stopcodon_array.reset_index(inplace=True, drop=True)
        threeprime_array = temp_df[temp_df.feature == '3\'-UTR']
        threeprime_array.reset_index(inplace=True, drop=True)

        fiveprime_intron_array = annotate_intron(fiveprime_array)
        cds_intron_array = annotate_intron(cds_array)
        threeprime_intron_array = annotate_intron(threeprime_array)
        transcript_with_introns_temp = pd.concat([transcript_array, fiveprime_intron_array, startcodon_array,
                                             cds_intron_array, stopcodon_array, threeprime_intron_array],
                                            ignore_index=True)
        transcript_with_introns = transcript_with_introns.append(transcript_with_introns_temp, ignore_index=True)
    transcript_with_introns.to_csv(args.intron_annotated_out, sep="\t", header=True, index=False)
    return transcript_with_introns


"""When the feature_array is empty line 98 throws and error"""
def annotate_intron(feature_array):
    intron_feature_array = pd.DataFrame(columns=col_names)
    seqname = feature_array.at[0, 'seqname']
    source = 'Custom_Intron_Script'
    feature = 'intron'
    score = '.'
    strand = feature_array.at[0, 'strand']
    frame = '.'
    attribute = feature_array.at[0, 'attribute']
    limit = (len(feature_array) - 1)
    if feature_array.empty:
        intron_feature_array = feature_array
    elif limit == 0:
        intron_feature_array = feature_array
    else:
        for n in range(0, (limit-1), 1):
            feature_entry = feature_array.loc[n]
            intron_start = int(feature_array.at[n, 'end'] + 1)
            intron_end = int(feature_array.at[n + 1, 'start'] - 1)
            intron_annotations = [seqname, source, feature, intron_start, intron_end, score, strand, frame,
                                  attribute]
            intron_entry = dict(zip(col_names, intron_annotations))
            intron_feature_array = intron_feature_array.append(feature_entry, ignore_index=True)
            intron_feature_array = intron_feature_array.append(intron_entry, ignore_index=True)
        last_feature_entry = feature_array.loc[limit]
        intron_feature_array = intron_feature_array.append(last_feature_entry, ignore_index=True)
    return intron_feature_array


def isnegeativegene(gene_chunck):
    return gene_chunck.at[1, 'strand'] in ['-', '- ']  # there is irradic spacing in the strand column
    # throughout the gtf files.


def reversechunck(gene_chunck):
    transcript_array = gene_chunck[gene_chunck.feature == 'transcript']
    features_array = gene_chunck[gene_chunck.feature != 'transcript']
    features_array = features_array.iloc[::-1]
    return transcript_array.append(features_array, ignore_index=True)


def meetsgenecriteria(featurelist):
    # return featurelist.count('5\'-UTR') >= 1 and featurelist.count('start_codon') == 1 and \
    #        featurelist.count('CDS') >= 1 and featurelist.count('stop_codon') == 1 and \
    #        featurelist.count('3\'-UTR') >= 1
    return featurelist.count('start_codon') == 1 and featurelist.count('CDS') >= 1 and \
           featurelist.count('stop_codon') == 1
    # The requirement for feature check of UTR was removed as mRNA can be leaderless in all domains of life.

with open('UTR-GTFs/test.gtf', "r") as file:
    completegenes = filter_genes(file)
    # Test if the filtering step was a success
    if completegenes.empty:
        print("Genes Filtered = Failure!")
    else:
        print("Genes Filtered = Success!" + "\n" + f"Output is saved to...")  # Need to insert a string for saved file

    annotated_genes = get_intron(completegenes)
    if annotated_genes.empty:
        print("Introns Annotated = Failure!")
    else:
        print("Introns Annotated = Success!" + "\n" + f"Output is saved to...")

