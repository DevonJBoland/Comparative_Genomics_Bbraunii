


from Bio import SeqIO
import pandas as pd
from Bio.SeqRecord import SeqRecord
import argparse as ap

parser = ap.ArgumentParser(description="Program to extract all BUSCO gene single copy hits")
parser.add_argument('proteindb', type=str, help="Protein file in fasta format (ex:genome.fasta)")
parser.add_argument('tsv', type=str, help="Hits file from BUSCO (initially name full_table.tsv")
parser.add_argument('species', type=str, help="Organism name for fasta header")
args = parser.parse_args()


def extract(a, b, c):
    busco_id = str(a)
    status = str(b)
    header = str(c)
    if status == "Complete":
        seq_list = []
        seq = record_dic[header].seq
        gene = seq[0:]
        record = SeqRecord(gene, args.species, "","")
        seq_list.append(record)
        SeqIO.write(seq_list, busco_id+"-"+args.species+"-tmp.fasta", "fasta")
    elif status == "Duplicated":
        log.write(busco_id+" was duplicated in "+args.species+"\n")
    elif status == "Missing":
        log.write(busco_id+" was missing from "+args.species+"\n")


with open(args.species+"-extract.log", "w+") as log:
    with open(args.tsv, "r") as hits:
        df = pd.read_csv(hits, delimiter='\t', skiprows=2)
        df_tuple = [tuple(row) for row in df.values]
        limit = len(df_tuple)
        record_dic = SeqIO.to_dict(SeqIO.parse(args.proteindb, "fasta"))
        for i in range(0, limit):
            extract(df_tuple[i][0], df_tuple[i][1], df_tuple[i][2])
