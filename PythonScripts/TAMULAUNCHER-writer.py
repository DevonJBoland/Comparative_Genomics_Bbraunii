#!/Users/devonboland/anaconda3/envs/Python-Bubbles/bin/python3

"""
Designed to print out 1,189 commands into a single file for use with tamualuncher on Grace
-Devon J. Boland, short simple script
"""

with open("aln.txt", "r") as names:
    output = open("commands.txt", "w+")
    for line in names:
        line = line.rstrip('\n')
        output.write(
            "iqtree2 -s /scratch/user/devonjboland/Phytozome/MAFFT/aln-"
            +line+" -m MFP -bb 1000 -T 5\n")
    output.close()
