import os
from collections import defaultdict
import copy
from re import A

## HOW TO RUN ME
# modify the required paths and then get the output table in the 
# QC folder inside 
# python /home/jmendietaes/programas/pipelines/ChIP/cluster/02e_ChipSummaryTable.py

# AIM
# To gather ONLY peak calling stats in the cell consensus scripts
#####################
## TO MODIFY BY THE USER

# basepath of all the outputs from the pipeline (bamfiles folder
# is in here)
basePath="/home/jmendietaes/data/2021/chip/allProcessed"

# outpath of the peak calling step in the pipeline
outpath=f"{basePath}/furtherAnalysis/08_projectRestart"

#######################
# If you didn't aligned and need to run this
# for i in  ../bamfiles/valid/01_firstTry/*bam; do 
#   ii=$(basename $i); ii=(${ii//\./ }); ii=${ii[0]}; 
#   touch "summary_"${ii}".txt"; 
#   echo -e "sample name\tfastq name\tread count\tmillions\nREAD \
#COUNTS \n\n${ii}\t${ii}\t0\t0"> "summary_"${ii}".txt"; done

# get peak summary files
summaryPath=f"{outpath}/peakCalling/MACS2/peaks/summary"
peakSumFiles = os.listdir(summaryPath)
# get peak consensus files
peakConsensusP = f"{outpath}/peakCalling/MACS2/consensusPeaks/bySameChip"

outTable = f"{outpath}/QC"
if not os.path.exists(outTable):
    os.makedirs(outTable)

# get information of all chips
allChip = [li.split('.')[0].split('_')[2].split('-')[0] for li in peakSumFiles if li.startswith('summary')]
allChip = sorted(list(set(allChip)))
allCell = [li.split('.')[0].split('_')[1].split('-')[0] for li in peakSumFiles if li.startswith('summary')]
allCell = sorted(list(set(allCell)))

fileAssoc = {}
for cell in allCell:
    fileAssoc[cell] = {}
for li in peakSumFiles:
    if li.startswith('summary'):
        chip = li.split('.')[0].split('_')[2].split('-')[0]
        cell = li.split('.')[0].split('_')[1].split('-')[0]
        fileAssoc[cell][chip] = li


consensusFiles = os.listdir(peakConsensusP)

# look for the information for the table
# define variables with the location of the interest headers
positionsPeak = [ "NUMBER OF NARROW PEAKS",
            "NUMBER OF BROAD PEAKS",
            ]
posIndexPeak = [None] * len(positionsPeak)

header = 'FileType\tsampleID\t'
header += 'narrowPeaks\tnarrowFRiP\tbroadPeaks\tbroadFRiP\t'
#header += 'epicPeaks\tepicFRiP\t'
header += 'narrowConsensusPeaks\tbroadConsensusPeaks\n'

allText = ""
for cell in allCell:
    for chip in allChip:
        if chip in fileAssoc[cell]:
            # now for each of the peak files
            fi2 = fileAssoc[cell][chip]

            content = 'peakfile\t'
            # id
            content += f"{'_'.join(fi2.split('.')[0].split('_')[1:])}\t"
            

            ## add rest of information (peak related)
            # open file
            with open(f'{summaryPath}/{fi2}', 'r') as f:
                data2=f.readlines()
            # get position of each location
            for npo, po in enumerate(positionsPeak):
                for nda, da in enumerate(data2):
                    if po in da:
                        posIndexPeak[npo] = nda

            if posIndexPeak[0] is not None:
                # header += f'NarrowPeaks\t'
                content += f'{data2[posIndexPeak[0]].split()[-1]}\t'
                # header += f'NarrowFRiP\t'
                content += f'{data2[posIndexPeak[0] + 4].split()[-1]}\t'
            else:
                content += "None\tNone\t"

            if posIndexPeak[1] is not None:
                # header += f'BroadPeaks\t'
                content += f'{data2[posIndexPeak[1]].split()[-1]}\t'
                # header += f'BroadFRiP\t'
                content += f'{data2[posIndexPeak[1] + 4].split()[-1]}\t'
            else:
                content += "None\tNone\t"
            
            #if posIndexPeak[2] is not None:
            #    # header += f'epicPeaks\t'
            #    content += f'{data2[posIndexPeak[2]].split()[-1]}\t'
            #    # header += f'epicFRiP\t'
            #    content += f'{data2[posIndexPeak[2] + 2].split()[-1]}\t'
            #else:
            #    content += "None\tNone\t"


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
with open(f"{outTable}/summaryTable_onlyChip.tsv", 'w') as fout:
    fout.write(header)
    fout.write(allText)
