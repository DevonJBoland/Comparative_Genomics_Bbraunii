#!/Users/devonboland/anaconda3/envs/Python-Bubbles/bin/python3

import pandas as pd

with open("Orthogroups.GeneCount.tsv", "r") as file:
    df = pd.read_csv(file, delimiter="\t")
    orthogroups = df.Orthogroup.to_list()
    species = list(df.columns.values[1:17])
    df.drop(labels=["Total"], axis=1, inplace=True)

# important_species = [3, 4, 5, 6, 8, 9]
# Create an empty list for Chlamy and braunii species
Chloso1602=[]
Chlre5_6=[]
Chlin1=[]
BobraA=[]
Botrbrau1=[]
BobraL=[]
for i in range(0,len(df.index)):
    if df.at[i, 'Chloso1602'] > 0:
        Chloso1602.append(df.at[i, "Orthogroup"])
    if df.at[i, 'Chlre5_6'] > 0:
        Chlre5_6.append(df.at[i, "Orthogroup"])
    if df.at[i, 'Chlin1'] > 0:
        Chlin1.append(df.at[i, "Orthogroup"])
    if df.at[i, 'BobraA'] > 0:
        BobraA.append(df.at[i, "Orthogroup"])
    if df.at[i, 'Botrbrau1'] > 0:
        Botrbrau1.append(df.at[i, "Orthogroup"])
    if df.at[i, 'BobraL'] > 0:
        BobraL.append(df.at[i, "Orthogroup"])


from matplotlib import pyplot as plt
from matplotlib_venn import venn3

Chloso1602 = set(Chloso1602)
Chlre5_6 = set(Chlre5_6)
Chlin1 = set(Chlin1)

venn3([Chloso1602, Chlre5_6, Chlin1], ('C. schcloesseri', 'C. reinhardtii', 'C. incerta')) # Assigns data, and circle names
plt.show() # show image

BobraA = set(BobraA)
Botrbrau1 = set(Botrbrau1)
BobraL = set(BobraL)

venn3([BobraA, Botrbrau1, BobraL], ('B. braunii A Race', 'B. braunii B Race', 'B. braunii L Race')) # Assigns data, and circle names
plt.show() # show image