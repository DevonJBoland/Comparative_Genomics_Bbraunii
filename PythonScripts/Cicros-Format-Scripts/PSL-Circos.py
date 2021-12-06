"""
Script to generate karyotupe and link file needed for input into Cricos, from PSL files
Created by Devon J. Boland
Created on: 08-31-2021
"""

import pandas as pd

with open("/Users/devonboland/Documents/Bbraunii-Omics/Cactus/halSynteny/Brace-Scaffolds-r.txt", "r") as Lkaryotype:
    Ldf = pd.read_csv(Lkaryotype, sep=' ', header=None, index_col=False)
    Lscaffolds = []
    for i in range(0, len(Ldf)):
        Lscaffolds.append(Ldf[2][i])

with open("/Users/devonboland/Documents/Bbraunii-Omics/Cactus/halSynteny/Arace-Scaffolds.txt", "r") as Akaryotype:
    Adf = pd.read_csv(Akaryotype, sep=' ', header=None, index_col=False)
    Ascaffolds = []
    for i in range(0, len(Adf)):
        Ascaffolds.append(Adf[2][i])

with open("/Users/devonboland/Documents/Bbraunii-Omics/Cactus/halSynteny/BvsA2.psl", "r") as LvsA:
    coordinates = pd.read_csv(LvsA, sep='\t', index_col=False, header=None)


chr1 = []
qs = []
qe = []
chr2 = []
ts = []
te = []
for i in range(0, len(coordinates)):
    for rowx in Ascaffolds:
        if coordinates[9][i] == rowx:
            for rowy in Lscaffolds:
                if coordinates[13][i] == rowy:
                    #for x in range(0, len(Adf)):
                        #if coordinates[9][i] == Adf[3][x]:
                    chr1.append(str(coordinates[9][i]))
                    qs.append(str(coordinates[10][i]))
                    qe.append(str(coordinates[11][i]))
                    chr2.append(str(coordinates[13][i]))
                    ts.append(str(coordinates[14][i]))
                    te.append(str(coordinates[15][i]))
colors = []
colors_tup = [tuple(row) for row in Adf.values]
for x in range(0, len(chr1)):
    for y in range(0, len(colors_tup)):
        if chr1[x] == colors_tup[y][2]:
            colors.append("color="+colors_tup[y][6])

output = pd.DataFrame(list(zip(chr1, qs, qe, chr2, ts, te, colors)), columns=["chr1", "chr1-start", "chr1-end"
                                                                              , "chr2", "chr2-start", "chr2-end",
                                                                              "color"])

with open("/Users/devonboland/Documents/Bbraunii-Omics/Cactus/halSynteny/BvsA-segdup.txt", "w") as log:
    tup = [tuple(row) for row in output.values]
    for x in range(0, len(tup)):
        log.write(tup[x][0]+" "+tup[x][1]+" "+tup[x][2]+" "+tup[x][3]+" "+tup[x][4]+" "+tup[x][5]+" "+tup[x][6]+"\n")



