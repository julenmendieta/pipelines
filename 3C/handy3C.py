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

# Read matrix from test format
def readMatrix(indir):
    '''Load matrix from txt to list of lists'''
    matrix_end=open(indir,'r')
    matrix = [[float(m) for m in ma.split()] for ma in matrix_end.read().split('\n')[:-1]]
    matrix_end.close()
    return matrix


def getLocusWithGenomeIntTadbit(hic_data, resol, locusCh, regionStart, regionEnd=False, 
                                wholeGenome=True, concatemersBin=False, bias=False):
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
        or list (bin coords)if want to give specific not continuous values
    :param False regionEnd: integer with end of interest bins range or
        nothing if values range provided in regionStart
    :param True wholeGenome: True if you want to get data from the whole genome 
        or False to go just for the chromosome
    :param False concatemersBin: If you want to normalise by frecuency introduce
        the dictionary (index in bin positions) with the number of times each bin 
        was found in a chimeric read
    :param False bias: TADbit bias dictionary
    :returns: A defaultdict with bins from the whole matrix as column
        and interactions floating or integer values
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

    # Get maximum bin in our data
    maxBin = 0
    for c in hic_data.section_pos:
        maxBin = max(maxBin, hic_data.section_pos[c][1])
    # Obtain number of cells in our data
    nCells = maxBin**2
    totalBins = maxBin # Dont need to ad +1 because last value does not exist, is to count all in range

    
    if wholeGenome == True:
        # in case we normalise given biases from OneD for example
        if bias != False:
            for cbin in binRange:
                # cbin moves us to the position where interactions happened in our loci/s of interest
                # ps will move us for all the other bins in the genome to look for interactions
                for ps in range(0, totalBins):
                    interVal = hic_data[cbin, ps] 
                    if interVal != 0 and cbin != ps:
                        interVal = interVal / bias[cbin] / bias[ps]
                        if not np.isnan(interVal):
                            interList[ps] += interVal
        # if we normalise by frequency 
        elif concatemersBin != False:
            for cbin in binRange:
                # cbin moves us to the position where interactions happened in our loci/s of interest
                # ps will move us for all the other bins in the genome to look for interactions
                for ps in range(0, totalBins):
                    if hic_data[cbin, ps] != 0 and cbin != ps:
                        divider = concatemersBin[cbin] + concatemersBin[ps] - hic_data[cbin, ps]
                        interVal = hic_data[cbin, ps] / float(divider)
                        interList[ps] += interVal
        # if we use raw data             
        else:
            for cbin in binRange:
                # cbin moves us to the position where interactions happened in our loci/s of interest
                # ps will move us for all the other bins in the genome to look for interactions
                for ps in range(0, totalBins):
                    interVal = hic_data[cbin, ps]
                    if interVal != 0 and cbin != ps:
                        interList[ps] += interVal
            
                
                
    else:
        # in case we normalise given biases from OneD for example
        if bias != False:
            for cbin in binRange:
                # cbin moves us to the position where interactions happened in our loci/s of interest
                # ps will move us for all the other bins in the genome to look for interactions
                for ps in range(chromBinBeg, chromBinEnd):
                    interVal = hic_data[cbin, ps]
                    if interVal != 0 and cbin != ps:
                        interVal = interVal / bias[cbin] / bias[ps]
                        if not np.isnan(interVal):
                            interList[ps] += interVal
                            
        # if we normalise by frequency 
        elif concatemersBin != False:
            for cbin in binRange:
                # cbin moves us to the position where interactions happened in our loci/s of interest
                # ps will move us for all the other bins in the chromosome to look for interactions
                for ps in range(chromBinBeg, chromBinEnd):
                    if hic_data[cbin, ps] != 0 and cbin != ps:
                        divider = concatemersBin[cbin] + concatemersBin[ps] - hic_data[cbin, ps]
                        interVal = hic_data[cbin, ps] / float(divider)
                        interList[ps] += interVal
        # if we use raw data     
        else:
            for cbin in binRange:
                # cbin moves us to the position where interactions happened in our loci/s of interest
                # ps will move us for all the other bins in the chromosome to look for interactions
                for ps in range(chromBinBeg, chromBinEnd):
                    interVal = hic_data[cbin, ps]
                    if interVal != 0 and cbin != ps:
                        interList[ps] += interVal
            
                    
                    
    # create variable for cases when foucs bin do not interact at all
    #this could happend in promoter capture HiC when a fragments spans at least 3 bins. In
    #this case it could be that the middle bind does not have any interaction at all because
    #the sequenced points are the borders of the restriction fragment
    delFocus = []
    for bi in binRange:
        if interList[bi] == 0 and hic_data[bi, bi] == 0:
            delFocus += [bi]
            
    if len(delFocus) != 0:
        message = 'Points removed from focus due to no interaction data '
        if wholeGenome == True:
            message += 'in the whole genome:\n %s' %delFocus
        else:
            message += 'in the whole chromosome:\n %s' %delFocus
        print message
    
    # update binRange
    binRange = list(set(binRange) - set(delFocus))
    
    
    return interList, binRange

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

def getNeighbourInteraction(data1, interList, locusCh, focus, dataType = 'matrix', 
                            percentaje = 95, wholeGenome=False, concatemersBin=False,
                            bias=False):
    '''
    :param data1: List of lists (dataType='matrix') or hic_data object (dataType='tadbit')
        with interaction data
    :param interList: A dictionary with bins from the whole matrix as column
        and interactions as values
    :param locusCh: Chromosome of focus point
    :param 'matrix' dataType: String indicating wether our data comes from a list of lists
        ('matrix') or a HiC_data object ('tadbit') from TADbit 
    :param 95 percentaje: Percentyle thresshold at wich to look for interactios
    :param False wholeGenome: Wether to look in the whole genome or not
    :param False concatemersBin: If you want to normalise by frecuency introduce
        the dictionary (index in bin positions) with the number of times each bin 
        was found in a chimeric read
    :param False bias: TADbit bias dictionary
    
    '''
    # To check if an interactions is significant, we should see if the neighbour
    #bins also interact with the focus region, and/or if some of the other TOP interacting
    #bins also interact with it
    # First get values distribution
    # with percentaje we set the percentyle were we will set our value thresshold

    limit = np.nanpercentile([inte for inte in interList.values() if inte != 0], percentaje)

    # Obtain all bins that have interactions above the thresshold limit with
    #our bin or bins of interest
    thresKeys = set()
    for k in interList.keys():
        if interList[k] > limit:
            thresKeys.add(k)

    # It could be that they where no interactions in the diagonal of our viewPoint
    #what could cause that wont be check in here, so if not present, we add it
    for fo in focus:
        thresKeys.add(fo)
    
    # Now we will just store the filtered interactions of the filtered bins, so
    #wont see any interaction with our oiriginal bin if it does not pass also
    #the filtering in the interacting bin
    # Iterate over each one of the filtered bins
    if dataType == 'matrix':
        highPeaks = {}
        # if freq normalise data
        if concatemersBin != False:
            for kPeak in sorted(list(thresKeys)):
                highPeaks[kPeak] = {}
                # Get new interaction limit for this bin
                kLimit = np.percentile([(inte / float(concatemersBin[kPeak] + concatemersBin[ni] - inte)) 
                                        for ni, inte in enumerate(data1[kPeak]) if inte != 0], percentaje) 
                newThres = thresKeys - set([kPeak])
                # check interactions of this peak with the others
                for nt in newThres:
                    # we will get with what kPeak interacts above the thresshold, and could be
                    #that the same thing interacts above the thresshold with kPeak, so later
                    #on the edge would be included twice, but is ok. It will help us to discern
                    #cases where binx has biny above thresshold but not the other way around
                    val1 = (data1[kPeak][nt] / float(concatemersBin[kPeak] + concatemersBin[nt] 
                                                     - data1[kPeak][nt]))
                    if val1 > kLimit:
                        highPeaks[kPeak][nt] = val1
                           
        elif bias != False:
            for kPeak in sorted(list(thresKeys)):
                highPeaks[kPeak] = {}
                # Get new interaction limit for this bin
                kLimit = np.nanpercentile([(inte / bias[ni] / bias[kPeak]) 
                                        for ni, inte in enumerate(data1[kPeak]) 
                                            if inte != 0], 
                                       percentaje) 
                newThres = thresKeys - set([kPeak])
                # check interactions of this peak with the others
                for nt in newThres:
                    # we will get with what kPeak interacts above the thresshold, and could be
                    #that the same thing interacts above the thresshold with kPeak, so later
                    #on the edge would be included twice, but is ok. It will help us to discern
                    #cases where binx has biny above thresshold but not the other way around
                    val1 = data1[kPeak][nt] / bias[kPeak] / bias[nt]
                    if val1 > kLimit:
                        highPeaks[kPeak][nt] = val1
        # if raw data
        else:
            for kPeak in sorted(list(thresKeys)):
                highPeaks[kPeak] = {}
                # Get new interaction limit for this bin
                kLimit = np.percentile([inte for inte in data1[kPeak] if inte != 0], percentaje) 
                newThres = thresKeys - set([kPeak])
                # check interactions of this peak with the others
                for nt in newThres:
                    # we will get with what kPeak interacts above the thresshold, and could be
                    #that the same thing interacts above the thresshold with kPeak, so later
                    #on the edge would be included twice, but is ok. It will help us to discern
                    #cases where binx has biny above thresshold but not the other way around
                    if data1[kPeak][nt] > kLimit:
                        highPeaks[kPeak][nt] = data1[kPeak][nt]
            
                    
    elif dataType == 'tadbit':
        
        # If we want to just check at chromosome level
        if wholeGenome == False:
            binBeg = data1.section_pos[locusCh][0]
            binEnd = data1.section_pos[locusCh][1]
        # If we want to go through the whole genome
        elif wholeGenome == True:
            # Get maximum bin in our data
            maxBin = 0
            for c in hic_data.section_pos:
                maxBin = max(maxBin, hic_data.section_pos[c][1])
            # Obtain number of cells in our data
            binBeg = 0
            binEnd = maxBin # Dont need to ad +1 because last value does not exist, is to count all in range


        highPeaks = {}
        # if freq norm data
        if concatemersBin != False:
            for kPeak in sorted(list(thresKeys)):
                highPeaks[kPeak] = {}
                # Get new interaction limit for this bin
                kLimit = np.percentile([(data1[kPeak, i] / float(concatemersBin[kPeak] + 
                                                            concatemersBin[i] - data1[kPeak, i]))  
                                        for i in range(binBeg, binEnd) 
                                        if data1[kPeak, i] != 0], percentaje)
                newThres = thresKeys - set([kPeak])
                # check interactions of this peak with the others
                for nt in newThres:
                    # we will get with what kPeak interacts above the thresshold, and could be
                    #that the same thing interacts above the thresshold with kPeak, so later
                    #on the edge would be included twice, but is ok. It will help us to discern
                    #cases where binx has biny above thresshold but not the other way around
                    val1 = (data1[kPeak,nt] / float(concatemersBin[kPeak] + concatemersBin[nt] 
                                                    - data1[kPeak,nt]))
                    if val1 > kLimit:
                        highPeaks[kPeak][nt] = val1
        
        # if normalised by TADbit
        elif bias != False:
            for kPeak in sorted(list(thresKeys)):
                highPeaks[kPeak] = {}
                # Get new interaction limit for this bin
                kLimit = np.nanpercentile([(data1[kPeak, i] / bias[kPeak] / bias[i])  
                                        for i in range(binBeg, binEnd) 
                                            if data1[kPeak, i] != 0], 
                                        percentaje)
                newThres = thresKeys - set([kPeak])
                # check interactions of this peak with the others
                for nt in newThres:
                    # we will get with what kPeak interacts above the thresshold, and could be
                    #that the same thing interacts above the thresshold with kPeak, so later
                    #on the edge would be included twice, but is ok. It will help us to discern
                    #cases where binx has biny above thresshold but not the other way around
                    val1 = data1[kPeak,nt] / bias[kPeak] / bias[nt]
                    if val1 > kLimit:
                        highPeaks[kPeak][nt] = val1
                        
        # if raw data
        else:
            for kPeak in sorted(list(thresKeys)):
                highPeaks[kPeak] = {}
                # Get new interaction limit for this bin
                kLimit = np.percentile([data1[kPeak, i]  for i in range(binBeg, binEnd) 
                                        if data1[kPeak, i] != 0], percentaje)
                newThres = thresKeys - set([kPeak])
                # check interactions of this peak with the others
                for nt in newThres:
                    # we will get with what kPeak interacts above the thresshold, and could be
                    #that the same thing interacts above the thresshold with kPeak, so later
                    #on the edge would be included twice, but is ok. It will help us to discern
                    #cases where binx has biny above thresshold but not the other way around
                    if data1[kPeak,nt] > kLimit:
                        highPeaks[kPeak][nt] = data1[kPeak,nt]

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

    return edgeList, highPeaks


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
