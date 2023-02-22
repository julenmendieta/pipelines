import os
from collections import defaultdict
import pandas as pd

# Path where the bam files are located. With IgG subfolder for controls and
# main subfolder for ChIP bams
inpath = "/scratch/julen/ChIP/allData/08_restartProject/bamFiles/all"
# Path and name to store output metadata file
outpath = "/scratch/julen/ChIP/allData/08_restartProject/ChromHMM/metadata/metada.tsv"

controls = {}
useFiles = []
files = os.listdir(f"{inpath}")
for fi in files:
    if fi.endswith('bam'):
        if 'IgG' in fi:
            controls[fi.split('_')[0]] = f"{fi}"
        else:
            useFiles += [f"{fi}"]


print(controls)
info = {'Cell':[], 'ChIP':[], 'ChipBam':[], 'controlBam':[]}
for fi in useFiles:
    cell, chip = fi.split('/')[-1].split('_', 1)
    chip = chip.split('_')[0].split('-')[0].split('.')[0]
    info['Cell'] += [cell]
    info['ChIP'] += [chip]
    info['ChipBam'] += [fi]
    info['controlBam'] += [controls[cell]]
    
df = pd.DataFrame.from_dict(info)

# Store metadata table 
df.to_csv(outpath, sep="\t", index=False, header=False)
print(outpath)