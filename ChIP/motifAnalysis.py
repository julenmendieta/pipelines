import os
import numpy as np
import subprocess

import pybedtools

import lxml.html as lh
import codecs

from scipy.stats import ks_2samp
from collections import defaultdict
import re
import pandas as pd
# plot related libraries
from matplotlib import pyplot as plt
import seaborn as sns
from matplotlib.cbook import get_sample_data


### set of functions to get motifs with HOMER
def getPeakCommand(species, coordsFile, outMotif, size, background, mtfLens,
                  threads=4, mknown=False):
    
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
    if mknown != False:
        mknown = f' -mknown {mknown}'
    else:
        mknown = ''

    motifCheck_cmd = 'findMotifsGenome.pl %s %s %s -size %s %s -len %s -p %s%s &> %s'
    motifCheck_cmd = motifCheck_cmd %(coordsFile, species,
                                     outMotif, size, background_,
                                     mtfLens, threads, mknown,
                                     f'{outMotif}/Homer.log')
    return motifCheck_cmd

def runMotif(motifCheck_cmd):
    # run job
    process = subprocess.Popen(motifCheck_cmd, stdout=subprocess.PIPE, shell=True)
    output, error = process.communicate()
    
def runCommand(cmd):
    # run job
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    output, error = process.communicate()
    
def motifsInAll(outpath, chip, df, size,
                seqlen, rnd, species, mtfLens, background,
               inParams, nCPU=1, mknown=False):
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
                                            outMotif, seqlen, background, mtfLens,
                                            threads=nCPU, mknown=mknown)

        # in the future when i have all the bakcground will do this in parallel
        # FOR THIS COORDSFILE WILL NEED TO HAVE A RANDOM NUMBER
        if nCPU != 1:
            process = subprocess.Popen(motifCheck_cmd, stdout=subprocess.PIPE, shell=True)
            output, error = process.communicate()

        inParams += [[motifCheck_cmd]]
    return inParams

        
def motifsBySubset(outpath, size, foldChanges, chip, cell1, cell2, df_new, posu, posd,
                                  seqlen, rnd, species, mtfLens, background, inParams,
                                  nCPU=1, mknown=False):
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
    outMotif2_noBg = f'{outpath}/by_FoldChange/{chip}/{cell2}-Yes-{foldChanges[1]}_{cell1}-No-{foldChanges[0]}_defaultBg'
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
                                            outMotif1, seqlen, coordsFile2, mtfLens,
                                            threads=nCPU, mknown=mknown)
        if nCPU != 1:
            process = subprocess.Popen(motifCheck_cmd1, stdout=subprocess.PIPE, shell=True)
            output, error = process.communicate()
        inParams += [[motifCheck_cmd1]]
    # cell 1 default background
    if not os.path.exists(f'{outMotif1_noBg}/homerResults.html'):
        print('Cell1 with default background')
        motifCheck_cmd1_noBg = getPeakCommand(species, coordsFile1, 
                                            outMotif1_noBg, seqlen, background, 
                                            mtfLens, threads=nCPU, mknown=mknown)
        if nCPU != 1:
            process = subprocess.Popen(motifCheck_cmd1_noBg, stdout=subprocess.PIPE, shell=True)
            output, error = process.communicate()
        inParams += [[motifCheck_cmd1_noBg]]

    # cell2
    if not os.path.exists(f'{outMotif2}/homerResults.html'):
        print('Cell2 with cell1 bg')
        motifCheck_cmd2 = getPeakCommand(species, coordsFile2, 
                                            outMotif2, seqlen, coordsFile1, 
                                            mtfLens, threads=nCPU, mknown=mknown)
        if nCPU != 1:
            process = subprocess.Popen(motifCheck_cmd2, stdout=subprocess.PIPE, shell=True)
            output, error = process.communicate()
        inParams += [[motifCheck_cmd2]]
    # cell 2 default background
    if not os.path.exists(f'{outMotif2_noBg}/homerResults.html'):
        print('Cell2 with default background')
        motifCheck_cmd2_noBg = getPeakCommand(species, coordsFile2, 
                                            outMotif2_noBg, seqlen, background, 
                                            mtfLens, threads=nCPU, mknown=mknown)
        if nCPU != 1:
            process = subprocess.Popen(motifCheck_cmd2_noBg, stdout=subprocess.PIPE, shell=True)
            output, error = process.communicate()
        inParams += [[motifCheck_cmd2_noBg]]

    return inParams


def motifsBySubset_atac(outpath, size, foldChanges, chip, cell1, cell2, df_new, 
                                    posu, posd,
                                  seqlen, rnd, species, mtfLens, background, inParams,
                                  atac, nCPU=1, mknown=False,
                                   minOverlap=10):
    '''
    :param atac: A dictionary with cell ID as key and ATAC peak coordinates dataframe 
        as value
    :param 10 minOverlap: Minimum overlap between ChiP and ATAC to mantain peak
    
    '''
    ### create temporal file with subset 1       
    outMotif1 = f'{outpath}/by_FoldChange/{chip}/{cell1}-Yes-{foldChanges[0]}_{cell2}-No-{foldChanges[1]}_cellBg'
    outMotif1_noBg = f'{outpath}/by_FoldChange/{chip}/{cell1}-Yes-{foldChanges[0]}_{cell2}-No-{foldChanges[1]}_defaultBg'
    if size == True:
        # remove intersecting coordinates
        b = pybedtools.BedTool(atac[cell1.split('_')[0]][[
                    'chr', 
                      'start', 
                      'end', 
                      'interval_id']].to_string(header=False, 
                                               index=False), from_string=True)
        a = pybedtools.BedTool(df_new.loc[posu, ['chr', 
                                                 'start', 
                                                 'end', 
                                                 'interval_id']].to_string(header=False, 
                                                                                          index=False), from_string=True)
        df_temp = a.intersect(b).to_dataframe()
        df_temp.rename(columns={'name':'interval_id',
                               'chrom':'chr'}, inplace=True)
        df_temp = df_temp[['interval_id', 'chr', 'start', 'end']]
        

        midpos = df_temp['start'] + ((df_temp['end'] - df_temp['start']) // 2)
        df_temp['start'] = midpos - (seqlen // 2)
        df_temp['end'] = midpos + (seqlen // 2)
        outMotif1 += f'/{seqlen}'
        outMotif1_noBg += f'/{seqlen}'

    else:
        # remove intersecting coordinates
        b = pybedtools.BedTool(atac[cell1.split('_')[0]][[
                    'chr', 
                      'start', 
                      'end', 
                      'interval_id']].to_string(header=False, 
                                               index=False), from_string=True)
        a = pybedtools.BedTool(df_new.loc[posu, ['chr', 
                                                 'start', 
                                                 'end', 
                                                 'interval_id']].to_string(header=False, 
                                                                                          index=False), from_string=True)
        df_temp = a.intersect(b).to_dataframe()
        df_temp.rename(columns={'name':'interval_id',
                               'chrom':'chr'}, inplace=True)
        df_temp = df_temp[['interval_id', 'chr', 'start', 'end']]
        
        outMotif1 += '/noLen'
        outMotif1_noBg += '/noLen'
        
    # Filter out too short peaks
    keepPos = (df_temp['end'] - df_temp['start']) >= minOverlap
    df_temp = df_temp[keepPos]
    df_temp['strand'] = '+'
    # create coordinate file 1
    coordsFile1 = f'{outpath}/temp_{cell1}_{rnd}.bed'
    df_temp.to_csv(path_or_buf=coordsFile1, sep='\t', 
                   index=False, header=False)



    ### create temporal file with subset 2    
    outMotif2 = f'{outpath}/by_FoldChange/{chip}/{cell2}-Yes-{foldChanges[1]}_{cell1}-No-{foldChanges[0]}_cellBg'
    outMotif2_noBg = f'{outpath}/by_FoldChange/{chip}/{cell2}-Yes-{foldChanges[1]}_{cell1}-No-{foldChanges[0]}_defaultBg'
    if size == True:
        # remove intersecting coordinates
        b = pybedtools.BedTool(atac[cell2.split('_')[0]][[
                    'chr', 
                      'start', 
                      'end', 
                      'interval_id']].to_string(header=False, 
                                               index=False), from_string=True)
        a = pybedtools.BedTool(df_new.loc[posd, ['chr', 
                                                 'start', 
                                                 'end', 
                                                 'interval_id']].to_string(header=False, 
                                                                                          index=False), from_string=True)
        df_temp = a.intersect(b).to_dataframe()
        df_temp.rename(columns={'name':'interval_id',
                               'chrom':'chr'}, inplace=True)
        df_temp = df_temp[['interval_id', 'chr', 'start', 'end']]
        

        midpos = df_temp['start'] + ((df_temp['end'] - df_temp['start']) // 2)
        df_temp['start'] = midpos - (seqlen // 2)
        df_temp['end'] = midpos + (seqlen // 2)
        outMotif1 += f'/{seqlen}'
        outMotif1_noBg += f'/{seqlen}'

    else:
        # remove intersecting coordinates
        b = pybedtools.BedTool(atac[cell2.split('_')[0]][[
                    'chr', 
                      'start', 
                      'end', 
                      'interval_id']].to_string(header=False, 
                                               index=False), from_string=True)
        a = pybedtools.BedTool(df_new.loc[posd, ['chr', 
                                                 'start', 
                                                 'end', 
                                                 'interval_id']].to_string(header=False, 
                                                                                          index=False), from_string=True)
        df_temp = a.intersect(b).to_dataframe()
        df_temp.rename(columns={'name':'interval_id',
                               'chrom':'chr'}, inplace=True)
        df_temp = df_temp[['interval_id', 'chr', 'start', 'end']]
        
        outMotif2 += '/noLen'
        outMotif2_noBg += '/noLen'
        
    # Filter out too short peaks
    keepPos = (df_temp['end'] - df_temp['start']) >= minOverlap
    df_temp = df_temp[keepPos]
    df_temp['strand'] = '+'
    coordsFile2 = f'{outpath}/temp_{cell2}_{rnd}.bed'
    df_temp.to_csv(path_or_buf=coordsFile2, sep='\t', 
                   index=False, header=False)

    ### Now we run them using each other as background
    # cell1
    if not os.path.exists(f'{outMotif1}/homerResults.html'):
        print('Cell1 with cell2 bg')
        motifCheck_cmd1 = getPeakCommand(species, coordsFile1, 
                                            outMotif1, seqlen, coordsFile2, mtfLens,
                                            threads=nCPU, mknown=mknown)
        if nCPU != 1:
            process = subprocess.Popen(motifCheck_cmd1, stdout=subprocess.PIPE, shell=True)
            output, error = process.communicate()
        inParams += [[motifCheck_cmd1]]
    # cell 1 default background
    if not os.path.exists(f'{outMotif1_noBg}/homerResults.html'):
        print('Cell1 with default background')
        motifCheck_cmd1_noBg = getPeakCommand(species, coordsFile1, 
                                            outMotif1_noBg, seqlen, background, 
                                            mtfLens, threads=nCPU, mknown=mknown)
        if nCPU != 1:
            process = subprocess.Popen(motifCheck_cmd1_noBg, stdout=subprocess.PIPE, shell=True)
            output, error = process.communicate()
        inParams += [[motifCheck_cmd1_noBg]]

    # cell2
    if not os.path.exists(f'{outMotif2}/homerResults.html'):
        print('Cell2 with cell1 bg')
        motifCheck_cmd2 = getPeakCommand(species, coordsFile2, 
                                            outMotif2, seqlen, coordsFile1, 
                                            mtfLens, threads=nCPU, mknown=mknown)
        if nCPU != 1:
            process = subprocess.Popen(motifCheck_cmd2, stdout=subprocess.PIPE, shell=True)
            output, error = process.communicate()
        inParams += [[motifCheck_cmd2]]
    # cell 2 default background
    if not os.path.exists(f'{outMotif2_noBg}/homerResults.html'):
        print('Cell2 with default background')
        motifCheck_cmd2_noBg = getPeakCommand(species, coordsFile2, 
                                            outMotif2_noBg, seqlen, background, 
                                            mtfLens, threads=nCPU, mknown=mknown)
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


def getHomerHtmlInfo(newResults, fullName=False,
                    returnBgPerce=False):
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
    backgPos = [nh for nh,h in enumerate(header) if h == '% of Background'][0]
    
    # and rest of table elements
    nameS = []
    pvalS = []
    logpvalS = []
    percInTargetS = []
    percInBackGS = []
    for t in tr_elements[1:]:
        if fullName:
            nameS += [t[namePos].text_content()]
        else:
            nameS += [t[namePos].text_content().split('/')[0]]
        pvalS += [float(t[pPos].text_content())]
        logpvalS += [float(t[logPPos].text_content())]
        percInTargetS += [float(t[percPos].text_content()[:-1])]
        percInBackGS += [float(t[backgPos].text_content()[:-1])]

    if returnBgPerce:
        return nameS, pvalS, logpvalS, percInTargetS, percInBackGS
    else:
        return nameS, pvalS, logpvalS, percInTargetS
    
## some functions to filter motifs by expression and name
def motifPassByExpre(cell, rnaTable, motif, motifToGene, minExp):
    '''
    Function that returns True or False indicating if protein of the motif is expressed or not
    :param cell: name of the cell or list with the name of the cells whose maximum expression we 
        will take
    :param rnaTable: ddataframe with expression values. Gene Id column must be named Gene
    :param motif: full motif ID
    :param motifToGene: dictionary indicatinf to which gene name the names in the motifs must
        be associated with. First we split by ( and check whole name, if no match we split
        by - and keep the first position, if no match again we split by : and keep all position.In 
        the last scenario all genes separated by : must be expressed. If after all the name is not 
        in motifToGene the program will crash
    :param minExp: minimum expression value to count a gene as expressed
    
    '''
    # get expression dictionary
    # cell can be a list of cells
    if isinstance(cell, list):
        expreDict = dict((k, v) for k,v in zip(rnaTable['Gene'], rnaTable[cell].max(axis=1)))
    else:
        expreDict = dict((k, v) for k,v in zip(rnaTable['Gene'], rnaTable[cell]))

    # get motif if criteria passed
    # split by (, and check, then by -, keep first, and the by : and check both names)
    motifn = motif.split('(')[0].lower()
    if motifn not in motifToGene:
        motifn = motifn.split('-')[0]
        if motifn not in motifToGene:
            motifn = motifn.split(':')
    keep = False
    if isinstance(motifn, list):
        keep = True
        for mo in motifn:
            # get added expression of genes inside
            expe = sum([expreDict.get(m, 0) for m in motifToGene[mo]])
            if expe >= minExp:
                keep = keep * True
            else:
                keep = False
    else:
        expe = sum([expreDict.get(m, 0) for m in motifToGene[motifn]])
        if expe >= minExp:
            keep = True
        else:
            keep = False
            
    return keep

#[g for g in rnaTable['Gene'] if 'atf' in g]
#rnaTable[rnaTable['Gene'] == 'Atf2']

def mergeMotifByFamilies(allMotifGenes, rnaTable=None, cell=False, minPval = 0.01,
                        returnFamilies=False):
    '''
    This function will check motifs that different by a number and merge them if their 
        expression value distribution is not different to the other members
        
    :param allMotifGenes: dictionary with motif names as key and a list of genes in value
    :param rnaTable: ddataframe with expression values. Gene Id column must be named Gene
    :param cell: name of the cell or list with the name of the cells whose maximum expression we 
        will take
    :param 0.01 minPval: minimum p-value for the KS test
    '''

    # now we will merge motifs by families
    allmotifs = allMotifGenes.keys()
    motifFamilies = defaultdict(list)
    for mo in allmotifs:
        family = re.sub("\\d[a-z]?$", "", mo)
        motifFamilies[family] += [mo]

    finalMotifGenes = {}
    familyMembers = {}
    # we go again, compare genes expression, and if same we merge
    for family in motifFamilies:
        if len( motifFamilies[family]) == 1:
            finalMotifGenes[motifFamilies[family][0]] = allMotifGenes[motifFamilies[family][0]]
            familyMembers[motifFamilies[family][0]] = [motifFamilies[family][0]]
        else:
            toMerge = []
            allGenes = []
            for mo in motifFamilies[family]:
                # get the expression of the genes in the other TF
                otherGenes = []
                for mt in [m for m in motifFamilies[family] if m != mo]:
                    otherGenes += list(allMotifGenes[mt])
                otherGenes = set(otherGenes)
                pos = rnaTable['Gene'].isin(otherGenes)
                geneExpOt = list(rnaTable.loc[pos, cell])

                # get expression in the specific TF
                genes = list(allMotifGenes[mo])
                pos = rnaTable['Gene'].isin(genes)
                geneExp = list(rnaTable.loc[pos, cell])

                if isinstance(rnaTable, pd.DataFrame):
                    # check by Kolmogorov-Smirnoff test if they are qual
                    nelems = len(geneExp)
                    nOther = len(otherGenes)
                    criticVal = 1.63*np.sqrt((nelems+nOther)/(nelems*nOther))

                    # first we try for greater than normal distribution centred in zero
                    kval_g, pval_g = ks_2samp(geneExpOt, 
                                    geneExp, alternative='two-sided')
                    if (kval_g > criticVal) and (pval_g <= minPval):
                        print(f'{mo} from family {family} has different expression distribution')
                        finalMotifGenes[mo] = set(allMotifGenes[mo])
                        familyMembers[mo] = [mo]
                    else:
                        toMerge += [mo]
                        allGenes += [genes]
                else:
                    toMerge += [mo]
                    allGenes += [genes]

            # now we merge the ones that have to be merged
            genes = set(sum([list(allMotifGenes[t]) for t in toMerge], []))
            if len(toMerge) == 1:
                finalMotifGenes[toMerge[0]] = genes
                familyMembers[toMerge[0]] = [toMerge[0]]
            else:
                finalMotifGenes[f'{family}_family'] = genes
                familyMembers[f'{family}_family'] = toMerge


    if returnFamilies:
        return finalMotifGenes, familyMembers
    else:
        return finalMotifGenes
        

def getMotifToGene():
    motifToGene = {'AMYB':['MYBL1'], 'AP-1':['Jun', 'Fos'], 'AP1':['Jun', 'Fos'], 'Atf1':['Atf1'], 
              'Atf2':['Atf2'], 'Atf2':['Atf2'], 'Atf4':['Atf4'], 'Atf7':['Atf7'],
              'BATF':['Batf'], 'BMYB':['Mybl2'], 'BORIS':['Ctcf'], 'Bach1':['Bach1'],
              'Bach2':['Bach2'], 'Bcl11a':['Bcl11a'], 'Bcl6':['Bcl6'], 'CDX4':['CDX4'],
              'CEBP':['Cebpe', 'Cebpa'], 'CRE':['Creb1'], 'CREB5':['Creb5'], 'CTCF':['Ctcf'], 
              'Cdx2':['Cdx2'], 'Chop':['Cebpe', 'Cebpa'], 'E2A':['Tcf3'], 'EHF':['EHF'],
              'ELF1':['Elf1'], 'ELF3':['Elf3'], 'Elf4':['Elf4'], 'ELF5':['Elf5'], 'ERG':['Erg'], 'ETS':['Ets1', 'Ets2'],
              'ETS1':['Ets1'], 'E-box':['noGene'], 'RUNX':['Runx1', 'Runx2', 'Runx3'], 'ETV1':['Etv1'],
              'ETV4':['Etv4'], 'EWS':['Ewsr1'], 'FLI1':['Fli1'], 'Elk1':['Elk1'], 'Elk4':['Elk4'], 'Etv2':['Etv2'],
              'FOXK1':['Foxk1'], 'FOXP1':['Foxp1'], 'Fos':['Fos'], 'Fosl2':['Fosl2'], 'Foxa2':['Hnf3b', 'Foxa2'],
              'Foxf1':['Foxf1', 'FKHL5', 'FREAC1'], 'Foxh1': ['Foxh1', 'Fast1'], 'Foxo3':['Foxo3'],
              'Fra1':['Fra1', 'Fosl1'], 'Fra2':['Fra2', 'Fosl2'], 'GABPA':['Gabpa'], 'GATA':['Gata1', 'Gata2',
                                                                                            'Gata3', 'Gata4', 'Gata5',
                                                                                            'Gata6'],
              'IR3':['Ir3'], 'IR4':['Ir4'], 'GATA3':['Gata3'], 'DR4':['tnfrsf10A', 'Trailr1'], 'SCL':['Tal1'],
              'GFY':['Gfy'], 'Gata1':['Gata1'], 'Gata2':['Gata2'], 'Gata4':['Gata4'], 'Gata6':['Gata6'],
              'HLF':['Hlf'], 'HOXB13':['Hoxb13'], 'Hoxa11':['Hoxa11'], 'Hoxa13':['Hoxa13'], 'Hoxc9':['Hoxc9'],
              'Hoxd10':['Hoxd10'], 'Hoxd11':['Hoxd11'], 'Hoxd13':['Hoxd13'], 'IRF1':['Irf1'], 'IRF2':['Irf2'],
              'IRF3':['Irf3'], 'IRF4':['Irf4'], 'IRF8':['Irf8'], 'ISRE': ['Irf1', 'Irf2', 'Irf3', 'Irf4', 'Irf5',
                                                                          'Irf6', 'Irf7', 'Irf8', 'Irf9'],
              'Jun':['Jun'], 'JunB':['Junb'], 'JunD':['Jund'], 'KLF1':['Klf1'], 'KLF14':['Klf14'], 'KLF3':['Klf3'],
              'KLF5':['Klf5'], 'KLF6':['Klf6'], 'MITF':['Mitf'], 'MYB':['Myb'], 'MafF':['Maff'], 'MafK':['MafK'],
              'NF-E2':['Nfe2'], 'NF1':['Nf1'], 'NFE2L2':['Nfe2l2'], 'NFIL3':['Nfil3'], 'Nrf2':['Nrf2'],
              'Oct6':['Oct6', 'Pou3f1'], 'PRDM1':['Prdm1'], 'Pu.1':['Sp1'], 'IRF': ['Irf1', 'Irf2', 'Irf3', 'Irf4', 'Irf5',
                                                                          'Irf6', 'Irf7', 'Irf8', 'Irf9'],
              'REST-NRSF':['Prickle1'], 'RUNX1':['Runx1'],'RUNX2':['Runx2'], 'Ronin':['Ronin', 'Thap11'],
              'SPDEF':['Spdef', 'Pdef'], 'STAT1':['Stat1'], 'Stat3':['Stat3'], 'STAT4':['Stat4'], 
               'STAT5':['Stat5'], 'STAT6':['Stat6'],
              'Smad4':['Smad4'], 'Sox6':['Sox6'], 'Sp2':['Sp2'], 'Sp5':['Sp5'], 'SpiB':['Spib'], 
               'Stat3+il21':['il21'], 'TEAD': ['Tead1', 'Tead2', 'Tead3', 'Tead4'], 'TEAD1':['Tead1'],
               'TEAD3':['Tead3'], 'TEAD4': ['Tead4'], 'TRPS1':['Trps1'], 'Unknown':['noGene'], 
               'YY1':['Yy1'], 'ZEB1':['Zeb1'], 'ZNF143|STAF':['Znf143', 'Staf'], 'c-Jun-CRE':['Jun'],
               'Atf3':['atf3'], 'ETS:E-box':['Ets1', 'Ets2']
              }
    # convert all to lowercase
    motifToGene2 = {}
    for k in motifToGene:
        motifToGene2[k.lower()] = [m.lower() for m in motifToGene[k]]
        
    motifToGene = motifToGene2
    return motifToGene

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
