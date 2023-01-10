import pandas as pd

# get a list with all the motif names and submit it to hgnc
# This will help us to associate a motif with its target genes
#grep "^>" /scratch/julen/ATAC/allData/02_firstATAC/TOBIAS/motifs/known_jaspar_t100.motifs | sed 's/>//g' > homer_vertebrates_known_motifs.txt
# for line in $(cat homer_vertebrates_known_motifs.txt); do 
#     mapLib=(${line//_/ }); 
#     mapLib=${mapLib[o]}; 
#     echo ${mapLib}
# done > homer_vertebrates_known_motifs_.txt
# mv homer_vertebrates_known_motifs_.txt homer_vertebrates_known_motifs.txt

hgncFile = '/home/julen/programas/annotations/motifTFlink/hgnc-symbol-check.csv'
motifsFile = '/scratch/julen/ATAC/allData/02_firstATAC/TOBIAS/motifs/known_jaspar_t100.motifs'

allNames = []
with open(motifsFile, 'r') as f:
    for line in f:
        if line.startswith('>'):
            allNames += [line.rstrip()[1:]]

df = pd.read_csv(hgncFile, sep=',', skiprows=1, )

dict2 = {'Input':[], 'name':[]}
noMatch = []
for tf in set(df['Input']):
    for motif in allNames:
        if motif.startswith(f"{tf}_"):
            df.loc[df['Input'] == tf, 'Approved symbol']
            genes = list( df.loc[df['Input'] == tf, 'Approved symbol'])
            if len(genes) == 0:
                noMatch += [motif]
            else:
                for ge in genes:
                    dict2['Input'] += [motif]
                    dict2['name'] += [ge]
df2 = pd.DataFrame.from_dict(dict2)
# get TF with no match
pos = df2['name'].isnull()

outp = '/home/julen/programas/annotations/motifTFlink/homer_vertebrates_known_Motif-hgnc_TOBIAS_valid.csv'
df2[pos == False].to_csv(outp, header=False, sep='\t', index=False)

outp = '/home/julen/programas/annotations/motifTFlink/homer_vertebrates_known_Motif-hgnc_TOBIAS_noMatch.csv'
df2[pos].to_csv(outp, header=False, sep='\t', index=False)