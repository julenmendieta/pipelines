#import argparse
import os
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
import argparse

def getTopModels(paths, topCor = {}, action='w', jobTime=False,
                outpath='./', dcut='', maxd1=''):
    '''
    paths list: List with the paths in which both the interaction matrix and the 
        optimisation output files are located OR two lists first the one containing
        the paths to a matrix and then the one containing the paths in which the 
        optimisation files associated with this amtrices are located (same order)
    
    
    '''
    if isinstance(paths[0], list):
        matrixfiles = []
        optimFiles = []
        for npath, p in enumerate(paths[0]):
            matrixfiles += [p]
            optimFiles += [paths[1][npath]]
    else:
        matrixfiles = paths
        optimFiles = paths

    if jobTime == False:
        jobTime = '02:00:00'  # HH:MM:SS
    orderedKeys = {}
    for npath, path in enumerate(optimFiles):
        print path
        cell = path.split('/')[-3]
        id1 = path.split('/')[-2]
        # Get parameters resume
        files = []
        for o in os.listdir(path):
            if o.startswith('Val') and o.endswith('txt'):
                files.append(path + o)

        if len(files) != 0:
            # Open them and check ranking for dcutoff and maxdist
            allResults = []
            for fi in files:
                with open(fi, 'r') as f:
                    for line in f:
                        if not line.startswith('#'):
                            tempi = [float(li) for li in line.split()[2:]]
                            # If we do have a correlation value
                            #if not np.isnan(tempi[-1]): 
                            allResults.append(tempi)
            # Sort by correlation value
            allResults = sorted(allResults, key=lambda x: x[4], reverse=True)

            # get top correlation
            print dcut, maxd1
            nres = 0
            while nres != None:
            #for a in allResults:
                a = allResults[nres]
                if not np.isnan(a[-1]):
                    if a[3] == dcut and a[0] == maxd1:
                        if not cell in topCor:
                            topCor[cell] = {}
                            orderedKeys[cell] = []
                        topCor[cell][id1] = a
                        orderedKeys[cell].append(id1)
                        nres = None
                    else:
                        nres += 1
                else:
                    nres += 1
                if nres >= len(allResults):
                    nres = None

    # prepare file to run modelling
    outmo = outpath + 'modellinParams.txt'
    # Now get file names for the models to be taken
    cmd = ''
    topsPath = []
    for cell in topCor:
        expPath = '/'.join(path.split('/')[:-3]) + '/%s/' %cell
        
        #for regi in sorted(topCor[cell].keys()):
        for regi in orderedKeys[cell]:
            inPath = expPath + regi + '/'
            
            cmd += matrixfiles[npath] + '\t'  # Matrix path
            cmd += '%s\t%s\t' %(topCor[cell][regi][1], topCor[cell][regi][2]) # low up
            cmd += '%s\t' %topCor[cell][regi][3]  # Dcut
            cmd += '%s\n' %topCor[cell][regi][0]  # maxDist

    # write file
    with open(outmo, action) as f:
        f.write(cmd)

    return topCor

def ploter(infodict, label1, label2='', id1='', minVal = 0.1,
            topd={}):
    fig = plt.figure(figsize=(20,5))
    keys = sorted(infodict.keys())
    allkeys = []
    allvalue = []
    for m in keys:
        for c in infodict[m]:
            if c >= minVal:
                allkeys.append(m)
                allvalue.append(c)

    if len(allvalue) != 0:
        pd_dat = pd.DataFrame(
                        {label1: allkeys,
                        'Correlation':allvalue})

        #sns.violinplot(x='Dcutoff value', y="Correlation",
        #                                data=pd_dat, order=keys,
        #                inner='stick')
        sns.stripplot(x=label1, y="Correlation",
                                        data=pd_dat, order=keys)

        # Get maximum correlation and draw line
        maxi = max(allvalue)
        plt.axhline(maxi, color='red')

        plt.legend()
        plt.title('Correlation found at each %s %s in %s' %(label1, label2, id1))
        plt.xlabel(label1)
        plt.ylabel('Correlation')
        #plt.show()
        plt.ylim(0.1, 1)
        
        # add dcutoff with top correlation to list
        maxd = [0,0]
        for dc in infodict:
            if max(infodict[dc]) > maxd[1]:
                maxd = [dc, max(infodict[dc])]
        if maxd[0] not in topd:
            topd[maxd[0]] = 1
        else:
            topd[maxd[0]] += 1

        # add line in top correlation
        posi = None
        for nal, al in enumerate(sorted(infodict.keys())):
            if al == maxd[0]:
                posi = nal
        plt.axvline(posi, color='blue', linestyle='--')

        plt.ioff()
    else:
        print 'No files with correlation in %s' %id1
##############################################################

def checkAll(outpath, inpaths, show_dcut=False, dcut=False, topModels=False,
                chunkSize=False):

    ###########################################
    # minimal value of correlation to be accepted
    minVal = 0.1

    # check if we defined number of models to build and chunk
    if isinstance(inpaths[0], list):
        matrixfiles = []
        optimFiles = []
        for npath, p in enumerate(inpaths[0]):
            matrixfiles += [p]
            optimFiles += [inpaths[1][npath]]
    else:
        matrixfiles = inpaths
        optimFiles = inpaths
   
    # Show distribution by distance cutoffs
    topd = {}
    if show_dcut == True:
        topCor = {}
        for path in optimFiles:
            
            id1 = '%s_%s' %(path.split('/')[-3], path.split('/')[-2])
            # Get parameters resume
            files = []
            for o in os.listdir(path):
                if o.startswith('Val') and o.endswith('txt'):
                    files.append(path + o)

            if len(files) != 0:
                # Open them and check ranking for dcutoff and maxdist
                allResults = []
                for fi in files:
                    with open(fi, 'r') as f:
                        for line in f:
                            if not line.startswith('#'):
                                tempi = [float(li) for li in line.split()[2:]]
                                # If we do have a correlation value
                                #if not np.isnan(tempi[-1]): 
                                allResults.append(tempi)
                # Sort by correlation value
                allResults = sorted(allResults, key=lambda x: x[4], reverse=True)

                # get top correlation
                for a in allResults:
                    if not np.isnan(a[-1]):
                        topCor[id1] = a

                ## get correlations for all dcuts
                dcuts = {}
                for a in allResults:
                    dcuts[a[3]] = []

                # Store correlation values associated with each metric
                for a in allResults:
                    dcuts[a[3]] += [np.nanmax([float('nan'), a[4]])]

                ## Plot dcutoff
                ploter(dcuts, label1='Distance cutoffs', minVal=minVal,
                        id1=id1, topd=topd)
            
            else:
                print 'No files in %s' %id1

        for dc in sorted(topd):
            print '%s\t%s' %(dc, topd[dc])
        plt.show()


    # show maxdist distributions of the model correlations at selected dcutoff
    # Next we should check in the selected dcutoff which is the best maxdist
    topd = {}
    if dcut != False:
        topCor = {}
        for path in optimFiles:
            id1 = '%s_%s' %(path.split('/')[-3], path.split('/')[-2])
            # Get parameters resume
            files = []
            for o in os.listdir(path):
                if o.startswith('Val') and o.endswith('txt'):
                    files.append(path + o)

            if len(files) != 0:
                # Open them and check ranking for dcutoff and maxdist
                allResults = []
                for fi in files:
                    with open(fi, 'r') as f:
                        for line in f:
                            if not line.startswith('#'):
                                allResults.append([float(li) for li in line.split()[2:]])
                # Sort by correlation value
                allResults = sorted(allResults, key=lambda x: x[4], reverse=True)

                # get top correlation
                for a in allResults:
                    if not np.isnan(a[-1]):
                        topCor[id1] = a

                ## get correlations for all maxdist
                maxd = {}
                for a in allResults:
                    maxd[a[0]] = []

                # Store correlation values associated with each metric
                for a in allResults:
                    if a[3] == dcut:
                        maxd[a[0]] += [np.nanmax([float('nan'), a[4]])]

                ## Plot maxdist
                ploter(maxd, label1='Maxdist', label2='with dcutoff %s'%dcut,
                                minVal=minVal, id1=id1, topd=topd)

            else:
                print 'No files in %s' %id1

        for md in sorted(topd):
            print '%s\t%s' %(md, topd[md])
        plt.show()



    # Code to retrieve top models for the given dcut and mxd1 combination and prepare scp command
    #to take then to my computer
    #dcut = 80.0
    #maxd1 = 140.0
    if topModels != False:
        topCor = {}
        topModels = topModels.split('_')  # dcut_maxdist_jobtime
        dcut = float(topModels[0])
        maxd1 = float(topModels[1])

        topCor = getTopModels(inpaths, topCor = {}, action='w',
                                outpath=outpath, 
                                dcut=dcut, maxd1=maxd1)
    return topCor
