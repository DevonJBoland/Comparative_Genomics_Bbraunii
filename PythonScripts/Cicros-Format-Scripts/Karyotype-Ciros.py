"""
Script to generate karyotupe file needed for input into Cricos, from genome files
Created by Devon J. Boland
Created on: 08-31-2021
"""

import Bio.SeqIO as IO

colors = ["vvlgrey", "vlgrey", "lgrey","grey", "dgrey", "vdgrey", "vvdgrey", "vvlred", "vlred", "lred", "red", "dred",
           "vdred","vvdred", "vvlgreen", "vlgreen", "lgreen", "green", "dgreen", "vdgreen", "vvdgreen", "vvlblue",
          "vlblue","lblue", "blue", "dblue", "vdblue", "vvdblue", "vvlpurple", "vlpurple", "lpurple", "purple",
          "dpurple","vdpurple", "vvdpurple", "vvlorange", "vlorange", "lorange", "orange", "dorange", "vdorange",
          "vvdorange","vvlyellow", "vlyellow", "lyellow","yellow", "dyellow", "vdyellow", "vvdyellow"]

"""Arace = IO.to_dict(IO.parse("/Users/devonboland/Documents/Bbraunii-Omics/Cactus/halSynteny/Genome/Arace-Nuclear-Genome-v1.fasta", "fasta"))
with open("/Users/devonboland/Documents/Bbraunii-Omics/Cactus/halSynteny/Arace-Scaffolds.txt", "w") as output:
    counts = 0
    for key in Arace.items():
        if len(key[1].seq) >= 1000000:
            output.write("chr - "+key[0]+" "+"A"+str(counts)+" "+"0 "+str(len(key[1].seq))+" "+str(colors[counts])+'\n')
            counts += 1"""

Brace = IO.to_dict(IO.parse("/Users/devonboland/Documents/Bbraunii-Omics/Cactus/halSynteny/Genome/Bbraunii_502_2.0.fa", "fasta"))
with open("/Users/devonboland/Documents/Bbraunii-Omics/Cactus/halSynteny/Brace-Scaffolds-all.txt", "w") as output:
    counts = 0
    for key in Brace.items():
        if len(key[1].seq) >= 1:
            output.write("chr - "+key[0]+" "+str(counts)+" "+"0 "+str(len(key[1].seq))+" "+"grey"+'\n')
            counts += 1

"""Lrace = IO.to_dict(IO.parse("/Users/devonboland/Documents/Bbraunii-Omics/Cactus/halSynteny/Genome/Lrace-Nuclear-Genome-v1.fasta", "fasta"))
with open("/Users/devonboland/Documents/Bbraunii-Omics/Cactus/halSynteny/Lrace-Scaffolds.txt", "w") as output:
    counts = 0
    for key in Lrace.items():
        if len(key[1].seq) >= 1000000:
            output.write("chr - " + key[0] + " " + "L" + str(counts) + " " + "0 " + str(len(key[1].seq)) + " " + str(
                colors[counts]) + '\n')
            counts += 1"""
