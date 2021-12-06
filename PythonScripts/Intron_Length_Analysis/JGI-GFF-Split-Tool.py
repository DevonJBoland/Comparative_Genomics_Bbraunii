"""
Script meant to search for user provided term: exon, CDS, gene (start-stop codon length)
Created by: Devon J. Boland
Created on: 08-23-2021
"""

import pandas as pd
import argparse as ap


parser = ap.ArgumentParser(description="Script to extract feature from gff files")
parser.add_argument("feature", type=str, help="Feature to extract (CDS, exon, gene)")
parser.add_argument("input", type=str, help="gff file to extract feature from")
parser.add_argument("output", type=str, help="Path to write output to")
args = parser.parse_args()

feature = args.feature
#print("What feature do you want to extract?:")
#feature = input()
if f'{feature}' == "CDS":
    print("Extracting CDSs from "+args.input+"...")
    with open(args.input, "r") as file:
    #with open("/Users/devonboland/Documents/Bbraunii-Omics/JGI-Phytozome/chlorophyta_download_gff
    # /Astpho2/Astpho2_GeneCatalog_genes_20111209.gff") as file:
        data = pd.read_csv(file, sep="\t", header=None, index_col=False)
        data_tup = [tuple(row) for row in data.values]
        limit = len(data_tup)
        with open(args.output, "w") as log:
        #with open("/Users/devonboland/Documents/Bbraunii-Omics/JGI-Phytozome/
        # chlorophyta_download_gff/Astpho2/test.gff", "w") as log:
            for i in range(0, limit):
                if str(data_tup[i][2]) == f"{feature}":
                    log.write(data_tup[i][0]+"\t"+data_tup[i][1]+"\t"+data_tup[i][2]+"\t"+str(data_tup[i][3])+"\t"+
                              str(data_tup[i][4])+"\t"+data_tup[i][5]+"\t"+data_tup[i][6]+"\t"+data_tup[i][7]+
                              "\t"+data_tup[i][8]+"\n")
elif f"{feature}" == "exon":
    print("Extracting exons from "+args.input+"...")
    with open(args.input, "r") as file:
    #with open("/Users/devonboland/Documents/Bbraunii-Omics/JGI-Phytozome/chlorophyta_download_gff
    # /Astpho2/Astpho2_GeneCatalog_genes_20111209.gff") as file:
        data = pd.read_csv(file, sep="\t", header=None, index_col=False)
        data_tup = [tuple(row) for row in data.values]
        limit = len(data_tup)
        with open(args.output, "w") as log:
        #with open("/Users/devonboland/Documents/Bbraunii-Omics/JGI-Phytozome/chlorophyta_download_gff
        # /Astpho2/test.gff", "w") as log:
            for i in range(0, limit):
                if str(data_tup[i][2]) == f"{feature}":
                    log.write(data_tup[i][0]+"\t"+data_tup[i][1]+"\t"+data_tup[i][2]+"\t"+str(data_tup[i][3])+
                              "\t"+str(data_tup[i][4])+"\t"+data_tup[i][5]+"\t"+data_tup[i][6]+"\t"+data_tup[i][7]+
                              "\t"+data_tup[i][8]+"\n")
elif f'{feature}' == "gene":
    print("Extracting genes from "+args.input+"...")
    with open(args.input, "r") as file:
    #with open("/Users/devonboland/Documents/Bbraunii-Omics/JGI-Phytozome/chlorophyta_download_gff/
    # Astpho2/Astpho2_GeneCatalog_genes_20111209.gff") as file:
        data = pd.read_csv(file, sep="\t", header=None, index_col=False, names=['seqname', 'source', 'feature',
                                                                                  'start', 'stop', 'score', 'strand',
                                                                                  'frame', 'attribute'])
        data_tup = [tuple(row) for row in data.values]
        limit = len(data_tup)
        check = data['feature'].str.contains('gene').any()
        if f'{check}' == "True":
            with open(args.output, "w") as log:
            #with open("/Users/devonboland/Documents/Bbraunii-Omics/JGI-Phytozome/chlorophyta_download_gff
            # /Astpho2/test.gff", "w") as log:
                for i in range(0, limit):
                    if str(data_tup[i][2]) == f"{feature}":
                        log.write(data_tup[i][0]+"\t"+data_tup[i][1]+"\t"+data_tup[i][2]+"\t"+str(data_tup[i][3])+
                                  "\t"+str(data_tup[i][4])+"\t"+data_tup[i][5]+"\t"+data_tup[i][6]+"\t"+data_tup[i][7]+
                                  "\t"+data_tup[i][8]+"\n")
        elif f'{check}' == "False":
            with open(args.output, "w") as log:
            #with open("/Users/devonboland/Documents/Bbraunii-Omics/JGI-Phytozome/chlorophyta_download_gff
            # Astpho2/test.gff","w") as log:
                # create empty list for each column in GFF file format
                seq_name = []
                source = []
                trait = [] # renamed as feature is already in use above
                start = []
                stop = []
                score = []
                strand = []
                frame = []
                attribute = []
                for i in range(0, limit):
                    if str(data_tup[i][2]) == "start_codon" or str(data_tup[i][2]) == "stop_codon":
                        seq_name.append(data_tup[i][0])
                        source.append(data_tup[i][1])
                        trait.append(data_tup[i][2])
                        start.append(data_tup[i][3])
                        stop.append(data_tup[i][4])
                        score.append(data_tup[i][5])
                        strand.append(data_tup[i][6])
                        frame.append(data_tup[i][7])
                        attribute.append(data_tup[i][8])

                # Create a data frame from lists above
                gene_df = pd.DataFrame(list(zip(seq_name, source, trait, start, stop, score, strand, frame, attribute)),
                                            columns=['seqname', 'source', 'feature', 'start', 'stop', 'score', 'strand',
                                                     'frame', 'attribute'])
                frequency = gene_df["attribute"].value_counts().to_dict()
                gene_freq_items = {str(k):int(v) for k, v in frequency.items()}
                sorted_frequency = sorted(gene_freq_items.items())
                sorted_frequency_tup = tuple(sorted_frequency)
                end = len(sorted_frequency_tup)
                paired = []
                for i in range(0, end):
                    if sorted_frequency_tup[i][1] == 2:
                        paired.append(sorted_frequency_tup[i][0])
                gene_df_tup =[tuple(row) for row in gene_df.values]
                end = len(gene_df_tup)
                # Recreate lists
                seq_name = []
                source = []
                trait = []  # renamed as feature is already in use above
                start = []
                stop = []
                score = []
                strand = []
                frame = []
                attribute = []
                for i in range(0, end):
                    if gene_df_tup[i][8] in paired:
                        seq_name.append(gene_df_tup[i][0])
                        source.append(gene_df_tup[i][1])
                        trait.append(gene_df_tup[i][2])
                        start.append(gene_df_tup[i][3])
                        stop.append(gene_df_tup[i][4])
                        score.append(gene_df_tup[i][5])
                        strand.append(gene_df_tup[i][6])
                        frame.append(gene_df_tup[i][7])
                        attribute.append(gene_df_tup[i][8])

                paired_df = pd.DataFrame(list(zip(seq_name, source, trait, start, stop, score, strand, frame, attribute)),
                                            columns=['seqname', 'source', 'feature', 'start', 'stop', 'score', 'strand',
                                                     'frame', 'attribute'])
                paired_df_tup = [tuple(row) for row in paired_df.values]
                end = len(paired_df_tup)
                for i in range (0, end, 2):
                    if i >= end:
                        break
                    start = str(paired_df_tup[i][3])
                    stop = str(paired_df_tup[i+1][4])
                    log.write(paired_df_tup[i][0]+'\t'+paired_df_tup[i][1]+'\t'+f'{feature}'+'\t'+start+'\t'
                              +stop+'\t'+paired_df_tup[i][5]+
                              '\t'+paired_df_tup[i][6]+'\t'+paired_df_tup[i][7]+'\t'+paired_df_tup[i][8]+'\n')
else:
    print("Input is not recognized, must use CDS, exon, or gene!")
    exit()
