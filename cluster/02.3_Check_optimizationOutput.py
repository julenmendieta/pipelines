#import argparse
import os
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
import argparse

def getTopModels(paths, topCor = {}, nmodels=1000, action='w', jobTime=False,
                outpath='./', chunkSize=1):

    orderedKeys = {}
    for path in paths:
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
    # Now get file names for the models to be taken
    cmd = ''
    topsPath = []
    for cell in topCor:
        expPath = '/'.join(path.split('/')[:-3]) + '/%s/' %cell
        for regi in orderedKeys[cell]:
            cmds = 'scp %s julen@172.16.54.11:/scratch/julen/PCHiC_assesment/modelos/%s/%s/'
            f = expPath + '%s/finalModel/%s_%s_*modelsAll' %(regi, cell, regi)
            #f = expPath + '%s/opt_LF%sUF%sC%sMdis%s_5000bp.models' %(
            #    regi, topCor[cell][regi][1], topCor[cell][regi][2], topCor[cell][regi][3],
            #    topCor[cell][regi][0])
            cmds2 = cmds %(f, cell, regi)
            topsPath.append(cmds2)

        # prepare file to run modelling
        outmo = outpath + 'modellinParams.txt'
        
        #for regi in sorted(topCor[cell].keys()):
        for regi in orderedKeys[cell]:
            inPath = expPath + regi + '/'
            # check for the matrix
            if 'FraserMin' in expPath:
                matPath = inPath + [o for o in os.listdir(inPath) 
                                    if o.startswith('MatrixNormMin')][0]
                if jobTime == False:
                    jobTime = '01:00:00'  # HH:MM:SS
            elif 'Fraser' in expPath:
                matPath = inPath + [o for o in os.listdir(inPath) 
                                    if o.startswith('MatrixNormFreqFiltrd')][0]
                if jobTime == False:
                    jobTime = '01:00:00'  # HH:MM:SS
            elif 'RAOvirt' in expPath:
                matPath = inPath + [o for o in os.listdir(inPath) 
                                    if o.startswith('Matrix')][0]
                if jobTime == False:
                    jobTime = '01:00:00'  # HH:MM:SS
            elif 'RAO' in expPath:
                matPath = inPath + [o for o in os.listdir(inPath) 
                                    if o.startswith('MatrixOneD')][0]
                if jobTime == False:
                    jobTime = '02:00:00'  # HH:MM:SS

            elif 'OneD' in expPath:
                matPath = inPath + [o for o in os.listdir(inPath) 
                                    if o.startswith('MatrixOneD')][0]
                if jobTime == False:
                    jobTime = '02:00:00'  # HH:MM:SS

            else:
                matPath = inPath + [o for o in os.listdir(inPath) 
                                    if o.startswith('Matrix')][0]
                if jobTime == False:
                    jobTime = '01:00:00'  # HH:MM:SS

            cmd += matPath + '\t'  # Matrix path
            cmd += '%s\t%s\t' %(topCor[cell][regi][1], topCor[cell][regi][2]) # low up
            cmd += '%s\t' %topCor[cell][regi][3]  # Dcut
            cmd += '%s\t' %topCor[cell][regi][0]  # maxDist
            cmd += '%s\t' %jobTime
            cmd += '%s\t' %nmodels
            cmd += '%s\n' %chunkSize

    # write file
    with open(outmo, action) as f:
        f.write(cmd)

    return topsPath, topCor

# ask for input
parser = argparse.ArgumentParser(description='')
parser.add_argument('-op','--outPath', help='Files_outPath',required=True)
parser.add_argument('-sdc','--dcutDistrib', help='show_dcutoff_corr_distrib',required=True)
parser.add_argument('-dc','--dcutoff', help='dcutoff_mxd_corr_distrib',required=True)
parser.add_argument('-tm','--topmodels', 
                        help='get_topModels_parameters',required=True)  # dcut_maxdist_jobtime
parser.add_argument('-nm','--nmodels',help='output_nmodels', required=False)
parser.add_argument('-chnk','--chunks',help='models_per_job', required=False)


# load input
args = parser.parse_args()
outpath=args.outPath
show_dcut=args.dcutDistrib
show_dcut = True if show_dcut == 'True' else False
dcut=args.dcutoff
dcut=False if dcut == 'False' else float(dcut)
topModels=args.topmodels
topModels=False if topModels == 'False' else topModels
nmodels=args.nmodels
chunkSize =args.chunks

###########################################
# minimal value of correlation to be accepted
minVal = 0.1

# check if we defined number of models to build and chunk
if nmodels != None:
    nmodels = int(nmodels)
else:
    nmodels = 1200
if chunkSize != None:
    chunkSize = int(chunkSize)
else:
    chunkSize = 1

#############################################

# we load the regions to check from the matrixList.txt file
matPaths= []
with open(outpath + 'matrixList.txt','r') as f:
	for line in f:
		matPaths += line.split()
paths = []
for ma in matPaths:
    paths.append('/'.join(ma.split('/')[:-1]) + '/')



#paths = [
#'/scratch/devel/jmendieta/deLaat/modelling/brain/reg1/',
#'/scratch/devel/jmendieta/deLaat/modelling/liver/reg1/',
#'/scratch/devel/jmendieta/deLaat/modelling/3T3/reg1/'
#]

# Show distribution by distance cutoffs
topd = {}
if show_dcut == True:
    topCor = {}
    for path in paths:
        
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
            fig = plt.figure(figsize=(20,5))
            keys = sorted(dcuts.keys())
            allkeys = []
            allvalue = []
            for m in keys:
                for c in dcuts[m]:
                    if c >= minVal:
                        allkeys.append(m)
                        allvalue.append(c)

            pd_dat = pd.DataFrame(
                            {'Dcutoff value': allkeys,
                            'Correlation':allvalue})

            #sns.violinplot(x='Dcutoff value', y="Correlation",
            #                                data=pd_dat, order=keys,
            #                inner='stick')
            sns.stripplot(x='Dcutoff value', y="Correlation",
                                            data=pd_dat, order=keys)

            # Get maximum correlation and draw line
            maxi = max(allvalue)
            plt.axhline(maxi, color='red')

            plt.legend()
            plt.title('Correlation found at each Distance cutoffs in %s' %id1)
            plt.xlabel('Distance cutoffs')
            plt.ylabel('Correlation')
            #plt.show()
            plt.ylim(0.1, 1)
            
            # add dcutoff with top correlation to list
            maxd = [0,0]
            for dc in dcuts:
                if max(dcuts[dc]) > maxd[1]:
                    maxd = [dc, max(dcuts[dc])]
            if maxd[0] not in topd:
                topd[maxd[0]] = 1
            else:
                topd[maxd[0]] += 1

            # add line in top correlation
            posi = None
            for nal, al in enumerate(sorted(dcuts.keys())):
                if al == maxd[0]:
                    posi = nal
            plt.axvline(posi, color='blue', linestyle='--')

            plt.ioff()
        else:
            print 'No files in %s' %id1

    for dc in sorted(topd):
        print '%s\t%s' %(dc, topd[dc])
    plt.show()


# show maxdist distributions of the model correlations at selected dcutoff
# Next we should check in the selected dcutoff which is the best maxdist
topd = {}
if dcut != False:
    for path in paths:
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

            ## get correlations for all maxdist
            maxd = {}
            for a in allResults:
                maxd[a[0]] = []

            # Store correlation values associated with each metric
            for a in allResults:
                if a[3] == dcut:
                    maxd[a[0]] += [np.nanmax([float('nan'), a[4]])]

            ## Plot maxdist
            fig = plt.figure(figsize=(20,5))
            keys = sorted(maxd.keys())
            allkeys = []
            allvalue = []
            for m in keys:
                for c in maxd[m]:
                    if c >= minVal:
                        allkeys.append(m)
                        allvalue.append(c)
            if len(allvalue) != 0:
                pd_dat = pd.DataFrame(
                                {'Dcutoff value': allkeys,
                                'Correlation':allvalue})

                #sns.violinplot(x='Dcutoff value', y="Correlation",
                #                                data=pd_dat, order=keys,
                #                inner='stick')

                sns.stripplot(x='Dcutoff value', y="Correlation",
                                            data=pd_dat, order=keys)

                # Get maximum correlation and draw line
                maxi = max(allvalue)
                plt.axhline(maxi, color='red')

                plt.legend()
                plt.title('Correlation found at each Maxdist with dcutoff %s in %s' %(dcut, id1))
                plt.xlabel('Maxdist')
                plt.ylabel('Correlation')
                #plt.show()
                plt.ylim(0.1, 1)

                # add dcutoff with top correlation to list
                maxd2 = [0,0]
                for md in maxd:
                    if len(maxd[md]) != 0:
                        if max(maxd[md]) > maxd2[1]:
                            maxd2 = [md, max(maxd[md])]
                if maxd2[0] not in topd:
                    topd[maxd2[0]] = 1
                else:
                    topd[maxd2[0]] += 1

                # add line in top correlation
                posi = None
                for nal, al in enumerate(sorted(maxd.keys())):
                    if al == maxd2[0]:
                        posi = nal
                plt.axvline(posi, color='blue', linestyle='--')

                plt.ioff()
            else:
                print 'No files with correlation in %s' %id1

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
    topModels = topModels.split('_')  # dcut_maxdist_jobtime
    dcut = float(topModels[0])
    maxd1 = float(topModels[1])
    jobTime = topModels[-1]

    # if different cell lines
    #difCel = False

    #topCor = {}
    # if difCel == True:
    #     alltPath = []
    #     for npath, path in enumerate(paths):
    #         if npath == 0:
    #             topsPath, topCor = getTopModels([path], topCor = topCor, nmodels=1000, 
    #                                 action='w', jobTime=jobTime, outpath=outpath)
    #         else:
    #             topsPath, topCor = getTopModels([path], topCor = topCor, nmodels=1000, 
    #                                 action='a', jobTime=jobTime, outpath=outpath)

    #         alltPath.append(topsPath)
    # else:
    #     topsPath, topCor = getTopModels(paths, topCor = topCor, nmodels=1000, action='w',
    #                                 jobTime=jobTime, outpath=outpath)
    

    topsPath, topCor = getTopModels(paths, topCor = {}, nmodels=nmodels, action='w',
                            jobTime=jobTime, outpath=outpath, chunkSize=chunkSize)

