import os
import numpy as np
import subprocess

import lxml.html as lh
import codecs

# plot related libraries
from matplotlib import pyplot as plt
import seaborn as sns
from matplotlib.cbook import get_sample_data


### set of functions to get motifs with HOMER
def getPeakCommand(species, coordsFile, outMotif, size, background, mtfLens,
                  threads=4):
    
    '''
    if background file defined you might need to check number of regions
        and if less use -h option
    
    '''
    if not os.path.exists(outMotif):
        os.makedirs(outMotif)

    if not size:
        size = 'given'
    if background:
        background_ = f'-bg {background}'
    else:
        background_ = ''
    if not mtfLens:
        mtfLens = '8,10,12'

    motifCheck_cmd = 'findMotifsGenome.pl %s %s %s -size %s %s -len %s -p %s &> %s'
    motifCheck_cmd = motifCheck_cmd %(coordsFile, species,
                                     outMotif, size, background_,
                                     mtfLens, threads,
                                     f'{outMotif}/Homer.log')
    return motifCheck_cmd

def runMotif(motifCheck_cmd):
    # run job
    process = subprocess.Popen(motifCheck_cmd2, stdout=subprocess.PIPE, shell=True)
    output, error = process.communicate()
    
    
def motifsInAll(outpath, chip, df, size,
                seqlen, rnd, species, mtfLens, background,
               inParams, nCPU=1):
    # create temporal file with all the peaks    
    outMotif = f'{outpath}/allPeakfile/{chip}/allPeak'
    if size == True:
        midpos = (df['end'] + df['start']) // 2
        df_temp = df.loc[:,['interval_id', 'chr', 'start', 'end']]
        df_temp['start'] = midpos - (seqlen // 2)
        df_temp['end'] = midpos + (seqlen // 2)
        outMotif += f'/{seqlen}'

    else:
        df_temp = df.loc[:,['interval_id', 'chr', 'start', 'end']]
        outMotif += '/noLen'
    df_temp['strand'] = '+'
    coordsFile = f'{outpath}/temp_{rnd}.bed'
    # contionue only if there are no results in there
    if not os.path.exists(f'{outMotif}/homerResults.html'):
        print('Running whole file')
        df_temp.to_csv(path_or_buf=coordsFile, sep='\t', 
                       index=False, header=False)
        motifCheck_cmd = getPeakCommand(species, coordsFile, 
                                            outMotif, seqlen, background, mtfLens)

        # in the future when i have all the bakcground will do this in parallel
        # FOR THIS COORDSFILE WILL NEED TO HAVE A RANDOM NUMBER
        if nCPU != 1:
            process = subprocess.Popen(motifCheck_cmd, stdout=subprocess.PIPE, shell=True)
            output, error = process.communicate()

        inParams += [[motifCheck_cmd]]
    return inParams

        
def motifsBySubset(outpath, size, foldChanges, chip, cell1, cell2, df_new, posu, posd,
                                  seqlen, rnd, species, mtfLens, background, inParams,
                                  nCPU=1):
    ### create temporal file with subset 1       
    outMotif1 = f'{outpath}/by_FoldChange/{chip}/{cell1}-Yes-{foldChanges[0]}_{cell2}-No-{foldChanges[1]}_cellBg'
    outMotif1_noBg = f'{outpath}/by_FoldChange/{chip}/{cell1}-Yes-{foldChanges[0]}_{cell2}-No-{foldChanges[1]}_defaultBg'
    if size == True:
        try:
            midpos = np.array([(
                    round(float(str(p).split(';')[0])) + s, 0) 
                        for p, s in zip(df_new.loc[posu, f'{cell1}.summit'], 
                                        df_new.loc[posu, 'start'])]) 
        except:
            nullPos = df_new.loc[:, f'{cell1}.summit'].isnull()
            print(f'{sum(posu & nullPos)} coordiantes had a peak in {cell1} but more signal in {cell2}')
            print(df_new.loc[posu & nullPos, ['chr', 'start', 'end', 'interval_id']])
            
            # change null values
            df_new = reasignMissingCols('summit', df_new, cell1, False, pos=(posu & nullPos))
            # get midpositions again
            midpos = np.array([(
                    round(float(str(p).split(';')[0])) + s, 0) 
                        for p, s in zip(df_new.loc[posu, f'{cell1}.summit'], 
                                        df_new.loc[posu, 'start'])]) 
            

        df_temp = df_new.loc[posu,['interval_id', 'chr', 'start', 'end']]
        df_temp['start'] = midpos - (seqlen // 2)
        df_temp['end'] = midpos + (seqlen // 2)
        outMotif1 += f'/{seqlen}'
        outMotif1_noBg += f'/{seqlen}'

    else:
        df_temp = df_new.loc[posu,['interval_id', 'chr', 'start', 'end']]
        outMotif1 += '/noLen'
        outMotif1_noBg += '/noLen'
    df_temp['strand'] = '+'
    # create coordinate file 1
    coordsFile1 = f'{outpath}/temp_{cell1}_{rnd}.bed'
    df_temp.to_csv(path_or_buf=coordsFile1, sep='\t', 
                   index=False, header=False)



    ### create temporal file with subset 2    
    outMotif2 = f'{outpath}/by_FoldChange/{chip}/{cell2}-Yes-{foldChanges[1]}_{cell1}-No-{foldChanges[0]}_cellBg'
    outMotif2_noBg = f'{outpath}/by_FoldChange/{chip}/{cell2}-Yes-{foldChanges[1]}_{cell1}-No-{foldChanges[0]}_noBg'
    if size == True:
        try:
            midpos = np.array([(
                    round(float(str(p).split(';')[0])) + s, 0) 
                        for p, s in zip(df_new.loc[posd, f'{cell2}.summit'], 
                                        df_new.loc[posd, 'start'])]) 
        except:
            nullPos = df_new.loc[:, f'{cell2}.summit'].isnull()
            print(f'{sum(posd & nullPos)} coordiantes had a peak in {cell2} but more signal in {cell1}')
            print(df_new.loc[posd & nullPos, ['chr', 'start', 'end', 'interval_id']])
            
            # change null values
            df_new = reasignMissingCols('summit', df_new, cell2, False, pos=(posd & nullPos))
            # get midpositions again
            midpos = np.array([(
                    round(float(str(p).split(';')[0])) + s, 0) 
                        for p, s in zip(df_new.loc[posd, f'{cell2}.summit'], 
                                        df_new.loc[posd, 'start'])]) 
            
        df_temp = df_new.loc[posd,['interval_id', 'chr', 'start', 'end']]
        df_temp['start'] = midpos - (seqlen // 2)
        df_temp['end'] = midpos + (seqlen // 2)
        outMotif2 += f'/{seqlen}'
        outMotif2_noBg += f'/{seqlen}'

    else:
        df_temp = df_new.loc[posd,['interval_id', 'chr', 'start', 'end']]
        outMotif2 += '/noLen'
        outMotif2_noBg += '/noLen'
    df_temp['strand'] = '+'
    coordsFile2 = f'{outpath}/temp_{cell2}_{rnd}.bed'
    df_temp.to_csv(path_or_buf=coordsFile2, sep='\t', 
                   index=False, header=False)

    ### Now we run them using each other as backgroun
    # cell1
    if not os.path.exists(f'{outMotif1}/homerResults.html'):
        print('Cell1 with cell2 bg')
        motifCheck_cmd1 = getPeakCommand(species, coordsFile1, 
                                            outMotif1, seqlen, coordsFile2, mtfLens)
        if nCPU != 1:
            process = subprocess.Popen(motifCheck_cmd1, stdout=subprocess.PIPE, shell=True)
            output, error = process.communicate()
        inParams += [[motifCheck_cmd1]]
    # cell 1 default background
    if not os.path.exists(f'{outMotif1_noBg}/homerResults.html'):
        print('Cell1 with default background')
        motifCheck_cmd1_noBg = getPeakCommand(species, coordsFile1, 
                                            outMotif1_noBg, seqlen, background, mtfLens)
        if nCPU != 1:
            process = subprocess.Popen(motifCheck_cmd1_noBg, stdout=subprocess.PIPE, shell=True)
            output, error = process.communicate()
        inParams += [[motifCheck_cmd1_noBg]]

    # cell2
    if not os.path.exists(f'{outMotif2}/homerResults.html'):
        print('Cell2 with cell1 bg')
        motifCheck_cmd2 = getPeakCommand(species, coordsFile2, 
                                            outMotif2, seqlen, coordsFile1, mtfLens)
        if nCPU != 1:
            process = subprocess.Popen(motifCheck_cmd2, stdout=subprocess.PIPE, shell=True)
            output, error = process.communicate()
        inParams += [[motifCheck_cmd2]]
    # cell 2 default background
    if not os.path.exists(f'{outMotif2_noBg}/homerResults.html'):
        print('Cell2 with default background')
        motifCheck_cmd2_noBg = getPeakCommand(species, coordsFile2, 
                                            outMotif2_noBg, seqlen, background, mtfLens)
        if nCPU != 1:
            process = subprocess.Popen(motifCheck_cmd2_noBg, stdout=subprocess.PIPE, shell=True)
            output, error = process.communicate()
        inParams += [[motifCheck_cmd2_noBg]]

    return inParams


# function to assign peak centre to summits where we dont have any data
def reasignMissingCols(whereCentre, df_new, cell1, cell2, pos=False):
    if not isinstance(pos, int):
        if whereCentre == 'summit':
            df_new.loc[pos, f'{cell1}.summit'] = (((df_new.loc[pos, 'start'] + 
                                                  df_new.loc[pos, 'end']) // 2) - 
                                                  df_new.loc[pos, 'start']) 
        elif whereCentre == 'midpos':
            df_new.loc[pos, f'{cell1}.summit'] = (((df_new.loc[pos, 'start'] + 
                                              df_new.loc[pos, 'end']) // 2) - 
                                              df_new.loc[pos, 'start']) 
        else:
            ERROR
        
    elif whereCentre == 'summit':
        pos1 = df_new.loc[:, f'{cell1}.summit'].isnull()
        df_new.loc[pos1, f'{cell1}.summit'] = (((df_new.loc[pos1, 'start'] + 
                                              df_new.loc[pos1, 'end']) // 2) - 
                                              df_new.loc[pos1, 'start']) 
        pos2 = df_new.loc[:, f'{cell2}.summit'].isnull()
        df_new.loc[pos2, f'{cell2}.summit'] = (((df_new.loc[pos2, 'start'] + 
                                              df_new.loc[pos2, 'end']) // 2) - 
                                              df_new.loc[pos2, 'start']) 
    elif whereCentre == 'midpos':
        ## cell 1
        # when peak no present we take consensus start
        pos11 = df_new.loc[:, f'{cell1}.start'].isnull()
        df_new.loc[pos11, f'{cell1}.summit'] = (((df_new.loc[pos11, 'start'] + 
                                              df_new.loc[pos11, 'end']) // 2) - 
                                              df_new.loc[pos11, 'start']) 
        # if peak present we take peak start
        pos12 = (df_new.loc[:, f'{cell1}.start'].isnull()) == False
        df_new.loc[pos12, f'{cell1}.summit'] = (((df_new.loc[pos12, f'{cell1}.start'] + 
                                              df_new.loc[pos12, f'{cell1}.end']) // 2) - 
                                              df_new.loc[pos12, f'{cell1}.start']) 

        ## cell 2
        # when peak no present we take consensus start
        pos2 = df_new.loc[:, f'{cell2}.start'].isnull()
        df_new.loc[pos2, f'{cell2}.summit'] = (((df_new.loc[pos2, 'start'] + 
                                              df_new.loc[pos2, 'end']) // 2) - 
                                              df_new.loc[pos2, 'start']) 
        # if peak present we take peak start
        pos22 = (df_new.loc[:, f'{cell2}.start'].isnull()) == False
        df_new.loc[pos22, f'{cell2}.summit'] = (((df_new.loc[pos22, f'{cell2}.start'] + 
                                              df_new.loc[pos22, f'{cell2}.end']) // 2) - 
                                              df_new.loc[pos22, f'{cell2}.start']) 
    else:
        ERROR

    return df_new


def getHomerHtmlInfo(newResults):
    f=codecs.open(newResults, 'r')
    doc = lh.fromstring(f.read())
    f.close()
    #Parse data that are stored between <tr>..</tr> of HTML
    tr_elements = doc.xpath('//tr')

    #Check the length of the first 12 rows
    if len(set([len(T) for T in tr_elements])) != 1:
        print('We got something appart than the table')
        ERROR

    # get header
    header=[]
    for t in tr_elements[0]:
        name=t.text_content()
        header.append(name)

    # locate the index of elements of interest
    namePos = [nh for nh,h in enumerate(header) if h == 'Best Match/Details'][0]
    pPos = [nh for nh,h in enumerate(header) if h == 'P-value'][0]
    logPPos = [nh for nh,h in enumerate(header) if h == 'log P-pvalue'][0]
    percPos = [nh for nh,h in enumerate(header) if h == '% of Targets'][0]

    # and rest of table elements
    nameS = []
    pvalS = []
    logpvalS = []
    percInTargetS = []
    for t in tr_elements[1:]:
        nameS += [t[namePos].text_content().split('/')[0]]
        pvalS += [float(t[pPos].text_content())]
        logpvalS += [float(t[logPPos].text_content())]
        percInTargetS += [float(t[percPos].text_content()[:-1])]

    return nameS, pvalS, logpvalS, percInTargetS
    

####################### PLOTS ########################
def plotMotifOccurrence(nameS, percInTargetS, pvalS=False, 
                        logpvalS=False, top=10, knownPath='',
                        title='', motifType='known', pdf=False):
    if pvalS == False and logpvalS == False:
        print('Nothing to plot')
    else:
        for npt, toPlot in enumerate([pvalS, logpvalS]):
            if toPlot != False:
                if npt == 0:
                    xlabel = 'Log 10 of P-value'
                else:
                    xlabel = '-Log e of P-value'
                    toPlot = [-1 * l for l in toPlot]
            
                fig = plt.figure(figsize=(20, 8))
                #ax = fig.add_subplot(2, 1, 1)
                #ax2 = fig.add_subplot(2, 2, 1)

                ax1 = plt.subplot2grid((top, 3), (0, 0), colspan=1, rowspan=top)
                _ = sns.barplot(y=nameS[:top], x=toPlot[:top], orient='h', 
                                   order=nameS[:top], ax=ax1)
                if npt == 0:
                    ax1.set_xscale('log', basex=10)
                ax1.set_xlabel(xlabel)

                # set axis limit
                ax1.set_xlim(min(toPlot[:top]) * 0.9, max(toPlot[:top]) * 1.01)

                ax2 = plt.subplot2grid((top, 3), (0, 1), colspan=1,  rowspan=top)
                line = sns.barplot(y=nameS[:top], x=percInTargetS[:top], orient='h', 
                                   order=nameS[:top], ax=ax2)
                ax2.set_xlim(0, 100)
                ax2.set_yticks([])
                ax2.set_xlabel('% of target sequences with motif')

                if title != '':
                    fig.suptitle(title)
                # add vertical line in times 2 concept of very significant value that
                # HOMER expects
                # Enrichment p-values reported by HOMER should be very very significant (i.e. << 1e-50)
                ax1.axvline(x=-np.log(1e-50) * 2, linestyle='--', color='red')

                # motif plots

                for i in range(0, top):
                    ax3 = plt.subplot2grid((top, 3), (i, 2), colspan=1)
                    motifLogo = f'{knownPath}/{motifType}{i+1}.logo.png'
                    try:
                        im = plt.imread(get_sample_data(motifLogo))
                        ax3.imshow(im)
                    except:
                        ax3.imshow([[]])
                    
                    ax3.axis('off')

                if pdf != False:
                    pdf.savefig(fig , bbox_inches='tight')
                plt.show()
