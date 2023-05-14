import pandas as pd
import os
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import sys

inpath = sys.argv[1]
files= os.listdir(inpath)
datFs = []
for fi in files:
    if fi.endswith("biotype_counts_mqc.tsv"):
        datFs += [pd.read_csv(f"{inpath}/{fi}", sep="\t", comment='#')]

df = datFs[0]
for i in range(len(datFs) - 1):
    df = pd.concat([df, datFs[i+1]])

# Get list of features and samples
cols = df.columns[1:]
samples = list(df["Sample"])

# Remove quote from sample names and percent from feature names
df['Sample'] = df['Sample'].str.replace(r"'", '')
df.rename(columns = dict((c, c[len("percent_"):]) for c in cols), 
        inplace=True)
df.index = df["Sample"]

# Plot and store
fig, ax = plt.subplots(1,1, figsize=(20, 5))
df.loc[sorted(df.index)].plot.barh(stacked=True, ax=ax)
plt.xlabel("% of Reads aligning to feature")
plt.legend(loc='center', bbox_to_anchor=(0.45, -0.4), ncol=7)
fig.suptitle("Gene type ")
plt.savefig(f"{inpath}/gene_types.pdf")
