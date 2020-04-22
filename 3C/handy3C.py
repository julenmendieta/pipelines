import pandas as pd
from collections import defaultdict
import numpy as np


#  Get a list of list from an interaction index file
def tableToMatrix(dirin, regionStart, regionEnd, resol=1, dirout=""):
    '''
    If positions from table already in bins from 0, leave resol
        as default, regionStart as 0 and regionEnd as the maximum
        bin
        Format of dirin would be: pos1  pos2    interaction
    '''
    end = False
    skip = 0
    # check for comments
    f=open(dirin, 'r')
    while end == False:
        line = f.readline()
        if line.startswith('#'):
            skip += 1
        else:
            end = True

    f=open(dirin, 'r')
    read_data = pd.read_csv(f,delim_whitespace=True,header=None,
                            skiprows=skip)

    # Tuple with the position of the matrix and value -> dict
    read_data['new_col'] = zip(read_data.ix[:,0], read_data.ix[:,1])
    dict_HiC = dict(zip([(k / resol, l / resol) for k, l in read_data.new_col], read_data[2]))
    f.close()

    len_chr = (regionEnd - regionStart) + resol # El ultimo bin tambien cuenta
    axis = (len_chr/resol)
    matrix = [[0 for _ in range(axis)]for _ in range(axis)]
    for i in range(axis):
        for j in range(axis):
            matrix[i][j] += dict_HiC.get((i + (regionStart/resol), j + (regionStart/resol)),0)
            if j != i:
                # Necessary to do a symetric part of the matrix
                matrix[j][i] += dict_HiC.get((i + (regionStart/resol), j + (regionStart/resol)),0)
    if dirout != "":
        matrix_end=open(dirout + "Matrix_%s_%s" %(regionStart, regionEnd),'w')
        matrix_end.write('\n'.join(['\t'.join([str(cell) for cell in line]) for line in matrix]) + '\n')
        matrix_end.close()
    else:
        return matrix


# Write matrix in text format
def writeMatrix(matrix, outdir):
    '''Write list of lists into matrix text file'''
    matrix_end=open(outdir,'w')
    matrix_end.write('\n'.join(['\t'.join([str(cell) for cell in line]) for line in matrix]) + '\n')
    matrix_end.close()

# Read matrix from text format
def readMatrix(indir):
    '''Load matrix from txt to list of lists'''
    matrix_end=open(indir,'r')
    matrix = [[float(m) for m in ma.split()] for ma in matrix_end.read().split('\n')[:-1]]
    matrix_end.close()
    return matrix


def normWayBias(interVal, cbin, ps, extraToNorm=None,
    positionAdjust=0):
    interVal = interVal / extraToNorm[cbin] / extraToNorm[ps]
    if not np.isnan(interVal):
        return interVal
    else:
        return 0
        
def normWayFreq(interVal, cbin, ps, extraToNorm=None,
    positionAdjust=0):
    divider = extraToNorm[cbin + positionAdjust] + \
                extraToNorm[ps + positionAdjust] - \
                interVal
    interVal = interVal / float(divider)
    return interVal
  
def normWayRaw(interVal, cbin, ps, extraToNorm=None,
    positionAdjust=0):
    return interVal

# from itertools
def chainxRange(*iterables):
    # chain('ABC', 'DEF') --> A B C D E F
    for it in iterables:
        for element in it:
            yield element

def getLocusWithGenomeIntTadbit(hic_data, resol, locusCh, regionStart, regionEnd=False, 
                                wholeGenome=True, concatemersBin=False, bias=False,
                                positionAdjust=0):
    '''
    Function to get all the interactions that bins from a given matrix have
        with a bin range of interest (from regionStart to regionEnd), or a
        list of positions (see bellow).
        Interactions are additive, so if we check bin 1 will have all
        the times it interacts with the bins of interest
    :param hic_data: A TADbit hic_data object
    :param resol: resolution of our experiment
    :param locusCh: Chromosome where our locus is located
    :param regionsStart: integer if starting point of a range (genomic coords),
        or list (bin coords)if want to give specific not continuous values. All
        positions must be from the same chromosome
    :param False regionEnd: integer with end of interest bins range or
        nothing if values range provided in regionStart
    :param True wholeGenome: True if you want to get data from the whole genome 
        or False to go just for the chromosome
    :param False concatemersBin: If you want to normalise by frecuency introduce
        the dictionary (index in bin positions) with the number of times each bin 
        was found in a chimeric read
    :param False bias: TADbit bias dictionary
    param 0 positionAdjust: In case the positions from concatemersBin are taking 
        into account bining from the whole genome, but we just load in hic_data
        one chromosome. Here concatemersBin will be added the 
        value in positionAdjust (which could also be negative) in order to 
        compensate this. In the example above, we would set positionAdjust
        to a negative number equal to the bin position in which the chromosome
        starts when using all data.

    :returns: One/two defaultdict with bins from the chromosome/chromosome+fullData
        as key and interaction floats or integers as values, depending on wholeGenome
        value False/True. It also resturns a list with the focus bins that had data
    '''
    chromBinBeg = hic_data.section_pos[locusCh][0]
    chromBinEnd = hic_data.section_pos[locusCh][1]
    # If we provide a range
    if isinstance(regionStart, int):
        binBeg = (chromBinBeg + (regionStart / resol))
        binEnd = (chromBinBeg + (regionEnd / resol))
        binRange = range(binBeg, binEnd + 1)
    # If instead provide the positions of interest in a list
    else:
        print 'Remember that focus location coordinates must be related to bins'
        binRange = regionStart
    

    #interList = {}
    interList = defaultdict(int)
    interChromList = defaultdict(int)

    # Get maximum bin in our data
    maxBin = 0
    for c in hic_data.section_pos:
        maxBin = max(maxBin, hic_data.section_pos[c][1])
    # Obtain number of cells in our data
    nCells = maxBin**2
    totalBins = maxBin # Dont need to ad +1 because last value does not exist, is to count all in range

    # set search start and end in case of whole genome search
    if wholeGenome == True:
        # is an iterator
        interRange = chainxRange(xrange(0, chromBinBeg), 
                                xrange(chromBinEnd+1, totalBins))
    else:
        interRange = []

    sstart = chromBinBeg
    send = chromBinEnd

    # taking into account type of comparison
    # in case we normalise given biases from OneD for example
    if bias != False:
        normWay = normWayBias
        extraToNorm = bias
     # if we normalise by frequency 
    elif concatemersBin != False:
        normWay = normWayFreq
        extraToNorm = concatemersBin
    # if we use raw data             
    else:
        normWay = normWayRaw
        extraToNorm = None
    
 
    for cbin in binRange:
        # cbin moves us to the position where interactions happened in our loci/s of interest
        # ps will move us for all the other bins in the genome to look for interactions
        for ps in xrange(sstart, send):
            interVal = hic_data[cbin, ps] 
            if interVal != 0 and cbin != ps:
                interVal = normWay(interVal, cbin, ps, 
                        extraToNorm=extraToNorm, 
                        positionAdjust=positionAdjust)
                interList[ps] += interVal

        for ps2 in interRange:
            interVal = hic_data[cbin, ps] 
            if interVal != 0 and cbin != ps:
                interVal = normWay(interVal, cbin, ps, 
                        extraToNorm=extraToNorm, 
                        positionAdjust=positionAdjust)
                interChromList[ps] += interVal

    # create variable for cases when focus bin do not interact at all
    #this could happend in promoter capture HiC when a fragments spans at least 3 bins. In
    #this case it could be that the middle bind does not have any interaction at all because
    #the sequenced points are the borders of the restriction fragment
    delFocus = []
    for bi in binRange:
        if interList[bi] == 0 and hic_data[bi, bi] == 0:
            delFocus += [bi]
    # in case we find them with interchrom
    if wholeGenome == True:
        reAd = []
        for d in delFocus:
            if interChromList[d] > 1:
                print 'Bin %s just have inter chromosomal interactions' %d
                reAd += [d]
        delFocus = list(set(delFocus) - set(reAd))
            
    if len(delFocus) != 0:
        message = 'Points removed from focus due to no interaction data, '
        message += '%s' %delFocus
        print message
    
    # update binRange
    binRange = list(set(binRange) - set(delFocus))
    
    if wholeGenome == True:
        return interList, binRange
    else:
        return interList, interChromList, binRange

def getLocusWithGenomeIntMatrix(matriz1, resol, regionStart, regionEnd):
    '''
    THIS FUNCTION HASNT BIN UPDATED AS THE ONE USING TADBIT DATA
    Function to get all the interactions that bins from a given matrix 
        (list of lists), WITH DIAGONAL DATA, have with a bin range of 
        interest (from regionStart to regionEnd), or a list of positions
        (see bellow).
        Interactions are additive, so if we check bin 1 will have all
        the times it interacts with the bins of interest
    :param matriz1: A matrix with the interactions data
    :param resol: resolution of our experiment
    :param regionsStart: integer if starting point of a range, or list
        if want to give specific not continuous values
    :param [] regionEnd: integer with end of interest bins range or
        nothing if values range provided in regionStart
    :returns: A defaultdict with bins from the whole matrix as column
        and interactions floating or integer values
    '''
    print 'THIS FUNCTION HASNT BIN UPDATED AS THE ONE USING TADBIT DATA'
    # If we provide a range
    if isinstance(regionStart, int):
        binBeg = regionStart / resol
        binEnd = regionEnd / resol
        binRange = range(binBeg, binEnd + 1)
    # If instead provide the positions of interest in a list
    else:
        binRange = regionStart

    #interList = {}
    interList = defaultdict(int)
    totalBins = len(matriz1) 
    
    for cbin in binRange:
        # cbin moves us to the position where interactions happened in our loci/s of interest
        # ps will move us for all the other bins in the chromosome to look for interactions
        for ps in range(0, totalBins):
            interVal = matriz1[cbin][ps]
            if interVal != 0:
                interList[ps] += interVal
                
    # create variable for cases when foucs bin do not interact at all
    #this could happend in promoter capture HiC when a fragments spans at least 3 bins. In
    #this case it could be that the middle bind does not have any interaction at all because
    #the sequenced points are the borders of the restriction fragment
    delFocus = []
    for bi in binRange:
        if interList[bi] == 0 and matriz1[bi][bi] == 0:
            delFocus += [bi]
    
    if len(delFocus) != 0:
        message = 'Points removed from focus due to no interaction data '
        message += 'in the whole matrix:\n %s' %delFocus
        print message
    
    # update binRange
    binRange = list(set(binRange) - set(delFocus))
    
                
    return interList, binRange


# Pair of functions to get top interactors
def computeRankPerc(inList, percentaje):
    return np.nanpercentile(inList, percentaje) 
                
def computeRankTop(inList, topInter):
    # in case list is sorter than said top
    topInter = min(topInter, len(inList)-1)
    return sorted(inList)[-topInter - 1]

def tadbitPeaks(data1, thresKeys, computeRank, 
                normWay, extraToNorm,
                positionAdjust, range1, topRank):
    highPeaks = {}
    # convert to list if generator
    if type(range1) != xrange:
        range1 = list(range1)
    for kPeak in sorted(list(thresKeys)):
        # have to create it constantly to avoid RAM usage
        highPeaks[kPeak] = {}
        # Get new interaction limit for this bin
        kLimit = computeRank([normWay(data1[kPeak, ni], kPeak, ni, 
                                        extraToNorm=extraToNorm,
                                        positionAdjust=positionAdjust)
                                for ni in range1  
                                    if data1[kPeak, ni] != 0], topRank)

        newThres = thresKeys - set([kPeak])
        # check interactions of this peak with the others
        for nt in newThres:
            # we will get with what kPeak interacts above the thresshold, and could be
            #that the same thing interacts above the thresshold with kPeak, so later
            #on the edge would be included twice, but is ok. It will help us to discern
            #cases where binx has biny above thresshold but not the other way around
            val1 = normWay(data1[kPeak,nt], kPeak, nt, extraToNorm=extraToNorm,
                            positionAdjust=positionAdjust)
            
            if val1 > kLimit:
                highPeaks[kPeak][nt] = val1
    
    return highPeaks

def getNeighbourInteraction(data1, interList, locusCh, focus, dataType = 'matrix', 
                            topRank = 95, wholeGenome=False, interChromList={},
                            concatemersBin=False,
                            bias=False, positionAdjust=0, topMeasure='percentyle'):
    '''
    :param data1: List of lists (dataType='matrix') or hic_data object (dataType='tadbit')
        with interaction data
    :param interList: A dictionary with bins from the whole matrix as column
        and interactions as values
    :param locusCh: Chromosome of focus point
    :param 'matrix' dataType: String indicating wether our data comes from a list of lists
        ('matrix') or a HiC_data object ('tadbit') from TADbit 
    :param 95 topRank: Percentyle thresshold at wich to look for interactios OR
        top value thresshold from which to look at interactions if topMeasure='top'
    :param False wholeGenome: Wether to look in the whole genome or not
    :param False concatemersBin: If you want to normalise by frecuency introduce
        the dictionary (index in bin positions) with the number of times each bin 
        was found in a chimeric read
    :param False bias: TADbit bias dictionary
    :param 0 positionAdjust: In case the positions from concatemersBin are taking 
        into account bining from the whole genome, but we just load in hic_data
        one chromosome. Here concatemersBin will be added the 
        value in positionAdjust (which could also be negative) in order to 
        compensate this. In the example above, we would set positionAdjust
        to a negative number equal to the bin position in which the chromosome
        starts when using all data.
    :param 'percentyle' topMeasure: Parameter to decide wether to use topRank to
        retrieve the percentyle of top interactors ('percentyle') or just the
        topRank top value in the list ('top'). Last one is not
        dependent of coverage
    
    '''

    # first define function to get maximum interactors
    if topMeasure == 'top':
        computeRank = computeRankTop
    elif topMeasure == 'percentyle':
        computeRank = computeRankPerc

    # taking into account type of comparison
    # in case we normalise given biases from OneD for example
    if bias != False:
        normWay = normWayBias
        extraToNorm = bias
    # if we normalise by frequency 
    elif concatemersBin != False:
        normWay = normWayFreq
        extraToNorm = concatemersBin
    # if we use raw data             
    else:
        normWay = normWayRaw
        extraToNorm = None

    # To check if an interactions is significant, we should see if the neighbour
    #bins also interact with the focus region, and/or if some of the other TOP interacting
    #bins also interact with it
    # First get values distribution
    # with topRank we set the percentyle  or top value were we will set our value thresshold
    
    limit = computeRank([inte for inte in interList.values() 
                                if inte != 0], topRank)
    # Obtain all bins that have interactions above the thresshold limit with
    #our bin or bins of interest
    thresKeys = set()
    for k in interList.keys():
        if interList[k] > limit:
            thresKeys.add(k)

    if wholeGenome == True:
        if interChromList == {}:
            print 'Missing interChromList to run whole genome'
            failll

        limit_interChr = computeRank([inte for inte in interChromList.values() 
                                        if inte != 0], topRank)
        thresKeysInter = set()
        for k in interChromList.keys():
            if interChromList[k] > limit_interChr:
                thresKeysInter.add(k)
        for fo in focus:
            thresKeysInter.add(fo)


    # It could be that they where no interactions in the diagonal of our viewPoint
    #what could cause that wont be check in here, so if not present, we add it
    for fo in focus:
        thresKeys.add(fo)
    
    # Now we will just store the filtered interactions of the filtered bins, so
    #wont see any interaction with our oiriginal bin if it does not pass also
    #the filtering in the interacting bin
    # Iterate over each one of the filtered bins
    if dataType == 'matrix':
        print 'Matrices are searched whole, no distinction for interChrom'
        highPeaks = {}
        
        for kPeak in sorted(list(thresKeys)):
            highPeaks[kPeak] = {}
            # Get new interaction limit for this bin
            kLimit = computeRank([normWay(inte, kPeak, ni, extraToNorm=extraToNorm,
                                            positionAdjust=positionAdjust)
                                    for ni, inte in enumerate(data1[kPeak]) 
                                        if inte != 0], topRank)

            newThres = thresKeys - set([kPeak])
            # check interactions of this peak with the others
            for nt in newThres:
                # we will get with what kPeak interacts above the thresshold, and could be
                #that the same thing interacts above the thresshold with kPeak, so later
                #on the edge would be included twice, but is ok. It will help us to discern
                #cases where binx has biny above thresshold but not the other way around
                val1 = normWay(data1[kPeak][nt], kPeak, nt, 
                                extraToNorm=extraToNorm,
                                positionAdjust=positionAdjust)
                
                if val1 > kLimit:
                    highPeaks[kPeak][nt] = val1    
                    
    elif dataType == 'tadbit':
        
        # at chromosome level
        chromBinBeg = data1.section_pos[locusCh][0]
        chromBinEnd = data1.section_pos[locusCh][1]
        intraRange = xrange(chromBinBeg, chromBinEnd)
        # If we want to go through the whole genome
        if wholeGenome != False:
            # Get maximum bin in our data
            maxBin = 0
            for c in data1.section_pos:
                maxBin = max(maxBin, data1.section_pos[c][1])
            # Obtain number of cells in our data
            binBeg = 0
            totalBins = maxBin # Dont need to ad +1 because last value does not exist, is to count all in range
            interRange = chainxRange(xrange(binBeg, chromBinBeg), 
                                xrange(chromBinEnd+1, totalBins))
            
            # get interChromosomal peaks
            highPeaksInter = tadbitPeaks(data1, thresKeysInter, computeRank, 
                                        normWay, extraToNorm,
                                        positionAdjust, interRange, topRank)

            ## Generate a list with all the interactions that passed the filtering
            edgeListInter = []
            for k in highPeaksInter.keys():
                if len(highPeaksInter[k]) != 0:
                    for k2 in highPeaksInter[k].keys():
                        edgeListInter.append((k, k2))

        # get intraChromosomal peaks
        highPeaks = tadbitPeaks(data1, thresKeys, computeRank, 
                                normWay, extraToNorm,
                                positionAdjust, intraRange, topRank)
        
    ## Generate a list with all the interactions that passed the filtering
    edgeList = []
    for k in highPeaks.keys():
        if len(highPeaks[k]) != 0:
            for k2 in highPeaks[k].keys():
                edgeList.append((k, k2))
                
    ## Generate list with all the interaction values
    #edgeVals = []
    #for k in highPeaks.keys():
    #    if len(highPeaks[k]) != 0:
    #        for k2 in highPeaks[k].keys():
    #            edgeVals.append(highPeaks[k][k2])
    #emini, emaxi =  min(edgeVals), max(edgeVals)
    if wholeGenome == False:
        return edgeList, highPeaks
    else:
        return edgeList, highPeaks, edgeListInter, highPeaksInter


def group_elements(element, group, visited, connections):
    # If this number has not been checked
    if visited[element]:
        return
    # Mark number as checked
    visited[element] = True
    # Add it to the groups list
    group.append(element)
    # Add the connected neigbourh numbers
    for another in connections[element]:
        group_elements(another, group, visited, connections)
    
    
    
def getInterGroups(edgeList, limitD = 10):
    '''
    Function to obtain the groups of conected bins
    :param edgeList: List with all the interactions between bins
    :param 10 limitD: Integer defining the maximum distance between two bins
        to set them as conected
    
    '''
    # Get the value distribution
    elements = set()
    for ed in edgeList:
        elements.add(ed[0])
        elements.add(ed[1])

    connections = {}
    visited = {}
    # For each number
    for element in elements:
        # Add the number to the conections list (empty)
        connections[element] = []
        # Add the number as non visited
        visited[element] = False
        # For each other numbers in list (itself included)
        for another in elements:
            # If the difference is smaller than the set limit
            #we add the other number as close interactor with 
            #first number (element)
            if abs(element - another) <= limitD:
                connections[element].append(another)

    groups = []
    for element in elements:
        if not visited[element]:
            group = []
            # Add the grouped elemnts to a list
            group_elements(element, group, visited, connections)
            groups.append(group)
    return sorted(groups)

def getGroupDegree(edgeList, groups):
    ## get degree of each group
    # First get appearances of each bin
    rept = []
    for ed in edgeList:
        rept += ed

    # Then turn appearances to degree
    degree = {}
    for response in rept:
        degree[response] = degree.setdefault(response, 0) + 1

    # Finally group degrees by group
    # for each group
    groupDegree = {}
    for n, g in enumerate(groups):
        # For each element inside
        groupSum = 0
        for gg in g:
            groupSum += degree[gg]
        # Store in a dictionary the sumatory of degrees and the group
        #members
        groupDegree[n] = {'degree' : groupSum, 'members' : sorted(g)}
    return groupDegree
