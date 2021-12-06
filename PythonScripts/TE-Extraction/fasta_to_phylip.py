from Bio import AlignIO

with open("consensi-algn.fa", 'r') as input:
    record = AlignIO.parse(input, "fasta")

    with open("test.phyl", "w") as output:
        AlignIO.write(record, output, "phylip")
