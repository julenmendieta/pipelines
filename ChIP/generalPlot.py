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




def TSSdistPlot(ax4, df_new, checkIndex):
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
            
            distanceToTSS = df_new.loc[checkIndex, 'Distance to TSS'].to_list()
            bottom = -sum(1 for d in distanceToTSS if d < 0)
            for di in TSSdists:
                presence = sum(1 for d in distanceToTSS if (di[0]*1000) <= d < (di[1] * 1000))
                if di[0] >= 0:
                    ax4.bar(['Distance to TSS'], [presence], 0.8, label='-'.join([str(d) for d in di]) + 'kb',
                          bottom=[bottom], color=TSSdists[di])
                else:
                    ax4.bar(['Distance to TSS'], [presence], 0.8, label='_nolegend_',
                          bottom=[bottom], color=TSSdists[di])


                bottom += presence

            ax4.set_ylabel('Count')
            #ax[2].set_title('Annotations presence in group')
            box = ax4.get_position()
            ax4.set_position([box.x0, box.y0, box.width * 0.8, box.height])

            # Put a legend to the right of the current axis
            ax4.legend(loc='upper center', bbox_to_anchor=(0.5, -0.08),
                  fancybox=False, shadow=False, ncol=3)

            ax4.hlines(y=[0], xmin=-0.5, xmax=0.5, colors='grey', alpha=0.5)
            
            
def motifPlot(groupIndexes, df2, motifsShort, axe, color_, topShowM=20):
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
                    20,
                  100)
    axe.set_xlabel(f'% of peaks with motif (from {len(groupIndexes)})')
    axe.vlines(x=[25, 50, 75], ymin=-0.5, ymax=topShowM, colors='grey', alpha=0.5)
    