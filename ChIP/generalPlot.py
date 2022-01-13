import pandas as pd
import os
from collections import defaultdict
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
import itertools
import matplotlib.backends.backend_pdf
sns.set(font_scale = 1.5)
sns.set_theme()




def TSSdistPlot(ax4, df_new, checkIndex, 
                bbox_to_anchor=(0.5, -0.08), showLegend=True,
                TSSrange='short'):
    '''
    Function to plot bar with distribution of distances to TSS
    :param 'short' TSSrange: plots only below or above 2.5 Kb as
        proximal or distal, respectively. Set to 'complete' for 
        a wider classification
    '''

    if TSSrange == 'complete':
        TSSdists = {(-float('inf'),-100):'#D795DC', 
                    (-100,-10):'#00C1B1', 
                    (-10,-5):'#76BB6D', 
                    (-5, -3):'#CCA662', 
                    (-3, -1):'#0083C2', 
                    (-1, 0):'#93CBE4',
                    (0,1):'#93CBE4', 
                    (1,3):'#0083C2', 
                    (3,5):'#CCA662', 
                    (5, 10):'#76BB6D', 
                    (10, 100):'#00C1B1', 
                    (100, float('inf')):'#D795DC'}
    else:
        TSSdists = {(-float('inf'),-2.5):'#CCA662', 
                    (-2.5, 0):'#93CBE4',
                    (0, 2.5):'#93CBE4', 
                    (2.5, float('inf')):'#CCA662'}
        labelDist = {(-float('inf'),-2.5): 'Distal',
                    (-2.5, 0): 'Proximal',
                     (0, 2.5): 'Proximal',
                     (2.5, float('inf')): 'Distal'}
    
    distanceToTSS = df_new.loc[checkIndex, 'Distance to TSS'].to_list()
    bottom = -sum(1 for d in distanceToTSS if d < 0)
    for di in TSSdists:
        presence = sum(1 for d in distanceToTSS if (di[0]*1000) <= d < (di[1] * 1000))
        if di[0] >= 0:
            if TSSrange == 'complete':
                ax4.bar(['Distance to TSS'], [presence], 0.8, 
                        label='-'.join([str(d) for d in di]) + 'kb',
                        bottom=[bottom], color=TSSdists[di], edgecolor = "none")
            else:
                ax4.bar(['Distance to TSS'], [presence], 0.8, 
                        label=labelDist[di],
                        bottom=[bottom], color=TSSdists[di], edgecolor = "none")
        else:
            ax4.bar(['Distance to TSS'], [presence], 0.8, label='_nolegend_',
                    bottom=[bottom], color=TSSdists[di], edgecolor = "none")


        bottom += presence

    ax4.set_ylabel('Count')
    # tight data limits
    ax4.set_ylim(-sum(np.array(distanceToTSS) < 0), sum(np.array(distanceToTSS) > 0))
    #ax[2].set_title('Annotations presence in group')
    box = ax4.get_position()
    ax4.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # Put a legend to the right of the current axis
    if showLegend:
        if TSSrange == 'complete':
            ax4.legend(loc='upper center', bbox_to_anchor=(0.5, -0.08),
                fancybox=False, shadow=False, ncol=3)
        else:
            ax4.legend(loc='upper center', bbox_to_anchor=(0.5, -0.08),
                fancybox=False, shadow=False, ncol=2)

    ax4.hlines(y=[0], xmin=-0.5, xmax=0.5, colors='grey', alpha=0.5)

# Plot annotation data    
def plotAnnot(df_start, ax, indexes, legend=True, allPosAnnt=None):
    if not allPosAnnt:
        allPosAnnt = {"3'":'#8da0cb',
                     "5'":'#fc8d62',
                     'Intergenic':'#b3b3b3',
                     'TTS':'#e78ac3',
                     'exon':'#a6d854',
                     'intron':'#66c2a5',
                     'non-coding':'#e5c494',
                     'promoter-TSS':'#ffd92f'}
    
    annot = df_start.loc[indexes]['Annotation'].str.split().str[0].sort_values().to_list()
    bottom = 0
    for an in allPosAnnt:
        if an in annot:
            presence = sum(1 for a in annot if a==an)
            ax.bar(['Annotation'], [presence], 0.8, label=an,
                 bottom=[bottom], color=allPosAnnt[an], edgecolor = "none")

            bottom += presence

    ax.set_ylabel('Count')
    # tight data limits
    ax.set_ylim(0, len(annot))
    #ax[3].set_title('Annotations presence in group')
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # Put a legend to the right of the current axis
    if legend:
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))


def motifPlot(groupIndexes, df2, motifsShort, axe, color_, topShowM=20,
                xlim=(20,100)):
    # get number of motifs in selected peaks
    noMotif = sum(df2.loc[groupIndexes][motifsShort].isnull().sum(axis=1) == len(motifsShort))
    motifSum = df2.loc[groupIndexes][motifsShort].notnull().sum()
    motifSum['NO MOTIF'] = noMotif
    motifSum.sort_values(ascending=True, inplace=True)

    # prepare plot of tendencies in the peak set outside groupIndexes
    motifSum_all = df2.loc[~df2.index.isin(groupIndexes)].copy()
    nonGroupIndexes = motifSum_all.index
    noMotif = sum(motifSum_all[motifsShort].isnull().sum(axis=1) == len(motifsShort))
    motifSum_all = motifSum_all[motifsShort].notnull().sum()
    motifSum_all['NO MOTIF'] = noMotif
    motifSum_all.sort_values(ascending=True, inplace=True)
    motifSumperce_all = motifSum_all / len(nonGroupIndexes) * 100
    motifSumperce_all.loc[motifSum.tail(topShowM).index].plot.barh(ax=axe, color='grey',
                                                                       width=0.7, alpha=0.7)

    # prepare plot of tendency in selection
    motifSumperce = motifSum / len(groupIndexes) * 100
    motifSumperce.tail(topShowM).plot.barh(ax=axe, color=color_, alpha=1,
                                          width=0.2)


    axe.set_title(f'Top {topShowM} motif occurrences. Selection vs rest (grey)')
    axe.set_xlim(#min(np.min(motifSumperce.tail(topShowM)) - (np.min(motifSumperce.tail(topShowM)) * 0.01), 70),
                    xlim[0],
                  xlim[1])
    axe.set_xlabel(f'% of peaks with motif (from {len(groupIndexes)})')
    axe.vlines(x=[25, 50, 75], ymin=-0.5, ymax=topShowM, colors='grey', alpha=0.5)
    