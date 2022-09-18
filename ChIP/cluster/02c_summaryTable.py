import os
from collections import defaultdict
import copy
from re import A

## HOW TO RUN ME
# modify the required paths and then get the output table in the 
# QC folder inside 
# python /home/jmendietaes/programas/PhD/ChIP/cluster/02c_summaryTable.py

#####################
## TO MODIFY BY THE USER
# path to the folder with all the sequencing QC
SequencedSum = "/home/jmendietaes/data/2021/chip/allProcessed/QC"

# basepath of all the outputs from the pipeline (bamfiles folder
# is in here)
basePath="/home/jmendietaes/data/2021/chip/allProcessed"

# outpath of the peak calling step in the pipeline
outpath=f"{basePath}/furtherAnalysis/subsampled_noIgG"

#######################


# get summary of sequencing of all dataset
seqSumFiles = [l for l in os.listdir(SequencedSum) if l.startswith('summary')]

# get peak summary files
summaryPath=f"{outpath}/peakCalling/MACS2/peaks/summary"
peakSumFiles = os.listdir(summaryPath)
# get peak consensus files
peakConsensusP = f"{outpath}/peakCalling/MACS2/consensusPeaks/bySameChip"

outTable = f"{outpath}/QC"
if not os.path.exists(outTable):
    os.makedirs(outTable)

## we need to know which files have been merged
# in the sequencing replicates
seqRep = os.listdir(f"{basePath}/bamfiles/valid/mergedSeqRep")
seqRep = [s.split('.sort')[0] for s in seqRep if s.endswith('bam')]
seqRepDict = {}
for chip in set([c.split('_')[1].split('-')[0] for c in seqRep]):
    seqRepDict[chip] = defaultdict(list)
    for s in seqRep:
        if f'_{chip}' in s:
            seqRepDict[chip][s.split('_')[0]] += [s]

# now real replicates
repli = os.listdir(f"{basePath}/bamfiles/valid/mergedReplicates")
repli = [s.split('.sort')[0] for s in repli if s.endswith('bam')]
repDict = {}
for chip in set([c.split('_')[1].split('-')[0] for c in repli]):
    repDict[chip] = defaultdict(list)
    for s in repli:
        if f'_{chip}' in s:
            repDict[chip][s.split('_')[0]] += [s]

# get information of all chips
allChip = [li.split('_')[2].split('-')[0] for li in seqSumFiles if li.startswith('summary')]
allChip = sorted(list(set(allChip)))

# Now we go through all final summaries and join them with 
# their sequencing summaries
# Special attention to subsampled cases
fileAssoc = {}
for chip in allChip:
    fileAssoc[chip] = {}
    allChIPCell = sorted(list(set([fi.split('_')[1] for fi in seqSumFiles if f'_{chip}' in fi])))
    for cell in allChIPCell:
        fileAssoc[chip][cell] = {'bamfiles':[], 'peakFiles':[]}
        for fi in seqSumFiles:
            if fi.startswith(f'summary_{cell}_{chip}'):
                id1 = '_'.join(fi.split('.')[0].split('_')[1:4])
                # look for the files associated with it
                for fi2 in peakSumFiles:
                    chip2 =  fi2.split('.')[0].split('_')
                    cell2 = chip2[1]
                    chip2 = chip2[2].split('-')[0]

                    # add connect it to the sequencing one
                    if (cell2 == cell) and (chip2 == chip):
                        fileAssoc[chip][cell]['bamfiles'] += [fi]
                        if not fi2 in fileAssoc[chip][cell]['peakFiles']:
                            fileAssoc[chip][cell]['peakFiles'] += [fi2]
# remove empty elements
fileAssoc2 = {}
for chip in fileAssoc:
    if len(fileAssoc[chip].keys()) != 0:
        for cell in fileAssoc[chip]:
            if len(fileAssoc[chip][cell]['bamfiles']) != 0:
                if not chip in fileAssoc2:
                    fileAssoc2[chip] = {}
                fileAssoc2[chip][cell] = fileAssoc[chip][cell]

fileAssoc = copy.copy(fileAssoc2)
consensusFiles = os.listdir(peakConsensusP)

# look for the information for the table
# define variables with the location of the interest headers
positions = [ "READ COUNTS \n",
            "SAMTOOLS FLAGSTAT - DUPLICATES \n",
            "SAMTOOLS FLAGSTAT - FINAL READS \n",
            ]
posIndex = [0] * len(positions)

positionsPeak = [ "NUMBER OF NARROW PEAKS",
            "NUMBER OF BROAD PEAKS",
            ]
posIndexPeak = [0] * len(positionsPeak)

header = 'fileType\tsampleID\ttotalReadPairs\ttotalMappedReadPairs\tmapedRatio\t'
header += 'duplicatedReadPairs\tdup%\tfinalReadPairs\tsubsampled\t'
header += 'narrowPeaks\tnarrowFRiP\tbroadPeaks\tbroadFRiP\t'
header += 'narrowConsensusPeaks\tbroadConsensusPeaks\n'

allText = ""
for chip in sorted(fileAssoc):
    for cell in sorted(fileAssoc[chip]):
        ## for different data merging scenarios
        # easiest one, everything was merged
        for fi in fileAssoc[chip][cell]['bamfiles']:
            with open(f'{SequencedSum}/{fi}', 'r') as f:
                data=f.readlines()

            # get position of each location
            for npo, po in enumerate(positions):
                for nda, da in enumerate(data):
                    if po == da:
                        posIndex[npo] = nda


            # get text to write  
            content = 'bamfile\t'
            # header = f'SampleID\t'
            content += f'{data[posIndex[0]+2].split()[0]}\t'

            # header += f'TotalReadPairs\t'
            content += f'{int(data[posIndex[0]+2].split()[2])*2}\t'

            # in few cases i dont have this data
            try:
                # header += f'TotalMappedReadPairs\t'
                # min 40Mill for TF and 80M for broad
                m = int(data[posIndex[1]+2].split()[0])
                content += f'{data[posIndex[1]+2].split()[0]}\t'

                # mapedRatio should be over 80% (in IgG around 60% ok)
                # The lower the dirtier our smaple is
                mr = int(data[posIndex[0]+2].split()[2])*2
                mr = round(m/float(mr), 2) * 100
                content += f'{mr}\t'

                # header += f'DuplicatedReadPairs\t'
                # should be < 20% the higher the less DNA content in sample
                content += f'{data[posIndex[1]+5].split()[0]}\t'

                # dup percentaje
                p = int(data[posIndex[1]+5].split()[0])
                content += f'{int(round((p/m) * 100,0))}\t'
            except:
                print(fi)
                content+='None\tNone\tNone\t'
        
            # header += f'FinalReadPairs\t'
            content += f'{data[posIndex[2]+2].split()[0]}\t'

            # header += f'Subsampled\t'
            content += f'False\t'

            allText += content + '\n'
            
        # now for each of the peak files
        for fi2 in fileAssoc[chip][cell]['peakFiles']:
            content = 'peakfile\t'
            # id
            content += f"{'_'.join(fi2.split('.')[0].split('_')[1:])}\t"
            # areas we cannot know without a lot of work
            content += '\t\t\t\t'
            if '-sub' in fi2.split('_')[2]:
                finalReads = fi2.split('.')[0].split('_')[2].split('sub')[-1]
                
                # header += f'FinalReadPairs\t'
                content += f'{finalReads}\t'

                # header += f'Subsampled\t'
                content += f'True\t'
            else:
                content += '\tFalse\t'

            ## add rest of information (peak related)
            # open file
            with open(f'{summaryPath}/{fi2}', 'r') as f:
                data2=f.readlines()
            # get position of each location
            for npo, po in enumerate(positionsPeak):
                for nda, da in enumerate(data2):
                    if po in da:
                        posIndexPeak[npo] = nda

            # header += f'NarrowPeaks\t'
            content += f'{data2[posIndexPeak[0]].split()[-1]}\t'
            # header += f'NarrowFRiP\t'
            content += f'{data2[posIndexPeak[0] + 4].split()[-1]}\t'

            # header += f'BroadPeaks\t'
            content += f'{data2[posIndexPeak[1]].split()[-1]}\t'
            # header += f'BroadFRiP\t'
            content += f'{data2[posIndexPeak[1] + 4].split()[-1]}\t'


            # Now consensus peaks info
            if (chip == 'IgG') or (chip == 'input'):
                nNarrowConsensus = 'None'
                nBroadConsensus = 'None'
            else:
                nNarrowConsensus = 0
                if f"{chip}_narrowPeak_consensusPeaks.bed" in os.listdir(peakConsensusP):
                    with open(f"{peakConsensusP}/{chip}_narrowPeak_consensusPeaks.bed", 'r') as ff:
                        for line in ff:
                            nNarrowConsensus += 1
                nBroadConsensus = 0
                if f"{chip}_broadPeak_consensusPeaks.bed" in os.listdir(peakConsensusP):
                    with open(f"{peakConsensusP}/{chip}_broadPeak_consensusPeaks.bed", 'r') as ff:
                        for line in ff:
                            nBroadConsensus += 1

            # header += f'narrowConsensusPeaks\t'
            content += f'{nNarrowConsensus}\t'

            # header += f'broadConsensusPeaks\t'
            content += f'{nBroadConsensus}\t'


            allText += content + '\n'

        allText += '\n'
    allText += '\n'

# Remove duplicates
at = []
for a in allText.split('\n'):
    if a != '':
        if a not in at:
            at += [a]
    else:
        at += [a]
allText = '\n'.join(at)

# write all
with open(f"{outTable}/summaryTable.tsv", 'w') as fout:
    fout.write(header)
    fout.write(allText)
