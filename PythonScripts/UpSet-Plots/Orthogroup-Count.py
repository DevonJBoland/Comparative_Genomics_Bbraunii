import pandas as pd


def strcheck(var, column):
    if str(var) != "nan":
        return names[column]
    else:
        return ""

names = ["BobraA","BobraL","BobraB","Chlorella_NC64A","Chlre5_6",
         "Coccomyxa_C163"]

#This block is meant to test the above function is working properly
#test_one = "this should be true"
#test_two = "nan"
#print(strcheck(test_one, 3))
#print(strcheck(test_two, 4))

workdir = "/Users/devonboland/Documents/Bbraunii-Omics/OrthoFinder-Analysis/Orthogroups/"

with open(workdir+"Orthogroups-bbraunii-chlorella-chlamy-cocco.txt", "r") as data:
    orthogroups = pd.read_csv(data, sep="\t", index_col=False)
orthogroups_tup = [tuple(row) for row in orthogroups.values]


with open (workdir+"upset-dataformat.tsv", "w") as log:
    for x in range(0, len(orthogroups_tup)):
        ortho = orthogroups_tup[x][0]
        a = strcheck(orthogroups_tup[x][1], 0)
        b = strcheck(orthogroups_tup[x][2], 1)
        c = strcheck(orthogroups_tup[x][3], 2)
        d = strcheck(orthogroups_tup[x][4], 3)
        e = strcheck(orthogroups_tup[x][5], 4)
        f = strcheck(orthogroups_tup[x][6], 5)
        #g = strcheck(orthogroups_tup[x][7], 6)
        #h = strcheck(orthogroups_tup[x][8], 7)
        #i = strcheck(orthogroups_tup[x][9], 8)
        #j = strcheck(orthogroups_tup[x][10], 9)
        #k = strcheck(orthogroups_tup[x][11], 10)
        #l = strcheck(orthogroups_tup[x][12], 11)
        #m = strcheck(orthogroups_tup[x][13], 12)
        #n = strcheck(orthogroups_tup[x][14], 13)
        #o = strcheck(orthogroups_tup[x][15], 14)
        #p = strcheck(orthogroups_tup[x][16], 15)
        log.write("\""+ortho+"\": "+"\""+a+","+b+","+c+","+d+","+e+","+f+"\""+","+"\n")
