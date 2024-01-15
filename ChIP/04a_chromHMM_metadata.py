import os
from collections import defaultdict
import pandas as pd

### Script to create metadata for ChromHMM
###
# Path where the bam files are located
#inpath = "/scratch/julen/ChIP/bamFiles/subsampled/all"
inpath="/scratch/julen/ChIP/bamFiles/subsampled/chromHMM"
#inpath="/scratch/julen/Arnau/allData/03_third_20231009/linkedBams/all"
# Path and name to store output metadata file
outpath = "/scratch/julen/ChIP/allData/08_restartProject/ChromHMM/metadata/metada.tsv"

controls = {}
useFiles = []
files = os.listdir(f"{inpath}")
for fi in files:
    if fi.endswith('bam'):
        if 'IgG' in fi:
        #if ('_H3-' in fi) or ('_H3_' in fi):
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