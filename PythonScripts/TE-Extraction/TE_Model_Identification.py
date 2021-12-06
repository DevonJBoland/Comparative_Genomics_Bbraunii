"""
Created by Devon J. Boland
Created on 06-28-2021

The purpose of this script is take the output from a blast (csv file) and extract the sequences
from the target alignment range, and output these alignments as a single fasta file.

The input for this script will be a list txt file containing the names for a set of files. These set of files will
be the sorted and dereplicated output of all the target alignments from the initial blastn.

Each file will be a csv file with three rows. The first row contains the fasta header for the location of the sequence
in the genome file. The second row will be the start of the alignment of the target respective to the query. The
third row will contain the end of the alignment of the target respective to the query.
"""
from Bio import SeqIO
import pandas as pd
from Bio.SeqRecord import SeqRecord
from contextlib import redirect_stdout
import argparse as ap

# Create Arguemnts Parser
parser = ap.ArgumentParser(description="Program to extract TE Model Alignments from a Genome file(input), using a "
                                       "post-formatted blastn output tsv file(hits) (outfmt6).The program will then "
                                       "extract the alignment plus x-bp(flank) of flanking sequence on either side. In the "
                                       "case that flanking sequence reaches either the 5' or 3' end of a sequence "
                                       "in the genome, then that end is the limit.")
parser.add_argument('input', type=str, help="Genome file in fasta format (ex:genome.fasta)")
parser.add_argument("hits", type=str, help="Hits file (ex:blast.outfmt6)")
parser.add_argument('flank', type=int, help='Length of flanking seqeunce to extract per side (ex:125)')
parser.add_argument("output", type=str, help="Output file (ex:extracted-tes.fasta)")
args = parser.parse_args()
seq_list=[]
def extract(a, b, c):
    # if statement to determine the actual start and stop regions of the alignmet in the + strand of the fasta file
    if b > c:
        d = c
        e = b
    elif c > b:
        d = b
        e = c
    header = str(a)
    start = int(d)
    end = int(e)
    seq = record_dic[header].seq
    if (start - args.flank) <= 0:
        start = 0
    if (start - args.flank) > 0:
        start = (start - args.flank)
    if (end + args.flank) >= len(seq):
        end = len(seq)
    if (end + args.flank) < len(seq):
        end = end + args.flank
    te_model = seq[start:end]
    record = SeqRecord(te_model, "TE_fragment_%i" % (i+1), "","")
    seq_list.append(record)
    with open("TE_Model_Extraction.log", "a+") as log:
        with redirect_stdout(log):
            print("The length of the TE_model_%i is" % (i+1), len(te_model), "while the length of the " \
                  "blastn alighnmet was", int(e - d))

with open(args.hits, "r") as targets:
    df = pd.read_csv(targets, delimiter='\t', header=None)
    blast_out_targets = [tuple(row) for row in df.values]
    record_dic = SeqIO.to_dict(SeqIO.parse(args.input, "fasta"))
    limit=len(blast_out_targets)
    for i in range(0,limit):
        extract(blast_out_targets[i][1], blast_out_targets[i][8], blast_out_targets[i][9])
SeqIO.write(seq_list, args.output, "fasta")
