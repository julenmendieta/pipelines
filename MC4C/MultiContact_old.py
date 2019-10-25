from collections import defaultdict
from pytadbit                     import HiC_data
import itertools

def factorial(n):return reduce(lambda x,y:x*y,[1]+range(1,n+1))


# Old function to normalise by frecuencies
def newNorm(hic_data, viewPoint, regionStartBin, regionEndBin, locusCh, tsvFile, 
            resol, method='', zeroPercent=False, mininter=0, keep=False, 
            returnNconcat=False):
    '''
    :param False keep: if you have a list of bin coordinates that should be the only 
        ones to be kept. e.j. [(1,3), (235,76), ...]
    '''
    dict_sec = hic_data.sections
    genome_seq = hic_data.chromosomes
    size = sum(genome_seq[crm] for crm in genome_seq)
    norm_data = HiC_data((), size, genome_seq, dict_sec, resolution=resol)
    
    
    sumStart = regionStartBin
    sumEnd = regionEndBin
    # Select till which point are we going to take into account interactions at
    #the time to normalise
    # Get iteration bins ranges
    regRange = set(range(regionStartBin, regionEndBin + 1))
    normRange = set(range(sumStart, sumEnd + 1))
    
    # Get concatemer appearances per bin
    if returnNconcat == False:
        concatemersBin = getConcatemersPerBin(hic_data, locusCh, regRange, tsvFile,
                                             resol, returnNconcat)
    else:
        concatemersBin, nConcat = getConcatemersPerBin(hic_data, locusCh, regRange, tsvFile,
                                         resol, returnNconcat)
    
    
    zeroSums = {}
    ## Filter data
    if zeroPercent != False:
        # For that, we turn it into zero
        # If badsThres is given in percentaje
        if zeroPercent <= 1:
            nbins = regionEndBin - regionStartBin + 1
            # Minimum number of non zero values
            badsThres = nbins * zeroPercent
        # If is given in minimum number of interactions (for one use 1.1)
        else:
            badsThres = zeroPercent
        bads = set()

        # First get 
        for bin1 in regRange:
                # If we are in the viewpoint just avoid diagonal
                zeroSums[bin1] = sum(i != bin1 and 
                                     not hic_data[bin1, i] and i not in viewPoint
                                     for i in normRange)
                if zeroSums[bin1] > badsThres:
                    bads.add(bin1)

                # get non zero cells
                suma = sum(i != bin1 and 
                            hic_data[bin1, i] and i not in viewPoint
                             for i in range(sumStart, sumEnd + 1))
                # If we have all zeros leave it as one for the enumerater, since
                #interaction value will alway be zero
                if suma == 0:
                    zeroSums[bin1] = 1

                else:
                    zeroSums[bin1] = ((zeroSums[bin1] + 1)/ float(suma))

        # Remove filtered columns
        regRange = regRange - bads
        normRange = normRange - bads
    
    # Now turn set to list again to recover order
    regRange = list(sorted(regRange))
    normRange = list(sorted(normRange))
    
    if zeroPercent == False:
        # get number of zeros if no filtering step
        for bin1 in regRange:
            # If we are in the viewpoint just avoid diagonal
            zeroSums[bin1] = sum(i != bin1 and 
                                 not hic_data[bin1, i] and i not in viewPoint
                                 for i in normRange)
            # If all points have data no adjustment required, so value = 1
            # pass

            suma = sum(i != bin1 and 
                            hic_data[bin1, i] and i not in viewPoint
                             for i in range(sumStart, sumEnd + 1))
            # If we have all zeros leave it as one for the enumerater, since
            #interaction value will alway be zero
            if suma == 0:
                zeroSums[bin1] = 1

            else:
                zeroSums[bin1] = ((zeroSums[bin1] + 1)/ float(suma))
                #zeroSums[bin1] = zeroSums[bin1]
                
    if method == 'totalConcat':
        if keep == False:
            for nbin1, bin1 in enumerate(regRange):
                for bin2 in regRange[nbin1:]:
                    # If diagonal or bellow sed interactions
                    if bin1 == bin2 or hic_data[bin1, bin2] <= mininter:
                        pass  # Leave it as zero

                    else:
                        # get divider
                        divider = nConcat
                        # if both are zero 
                        norm_data[bin1, bin2] = hic_data[bin1, bin2] / float(divider)
                        norm_data[bin2, bin1] = norm_data[bin1, bin2]
                    
        else:
            for ke in keep:
                # If diagonal or bellow sed interactions
                if (ke[0] == ke[1] or hic_data[ke[0], ke[1]] <= mininter or
                    (ke[0] not in regRange or ke[1] not in regRange)):
                    pass  # Leave it as zero
                else:
                    # get divider
                    divider = nConcat
                    # if both are zero 
                    norm_data[ke[0], ke[1]] = hic_data[ke[0], ke[1]] / float(divider)
                    norm_data[ke[1], ke[0]] = norm_data[ke[0], ke[1]]
        

    if method == 'zeroAdjust':
        for nbin1, bin1 in enumerate(regRange):
            for bin2 in regRange[nbin1:]:
                if bin1 == bin2:
                    norm_data[bin1, bin2] = 0
                elif hic_data[bin1, bin2] <= mininter:
                    norm_data[bin1, bin2] = 0
                elif bin1 in viewPoint or bin2 in viewPoint:
                    if concatemersBin[bin1] == 0:
                        if concatemersBin[bin2] == 0:
                            divider = 1
                            forZeros = 1
                        else:
                            divider = concatemersBin[bin2]
                            forZeros = zeroSums[bin2]
                    elif concatemersBin[bin2] == 0:
                        divider = concatemersBin[bin1]
                        forZeros = zeroSums[bin1]
                    else:
                        divider = concatemersBin[bin1] + concatemersBin[bin2]
                        forZeros = zeroSums[bin1] + zeroSums[bin2]
                    norm_data[bin1, bin2] = (hic_data[bin1, bin2] * divider**(1/2)) / float(forZeros)
                    norm_data[bin2, bin1] = norm_data[bin1, bin2]
                #elif bin2 in viewPoint and not bin1 in viewPoint:
                #    divider = float(concatemersBin[bin2]**2)
                #    forZeros = zeroSums[bin2]
                #    norm_data[bin1, bin2] = hic_data[bin1, bin2] / float(divider*forZeros)
                #    norm_data[bin2, bin1] = norm_data[bin1, bin2]
                else:
                    # get divider
                    if concatemersBin[bin1] == 0:
                        if concatemersBin[bin2] == 0:
                            divider = 1
                            forZeros = 1
                        else:
                            divider = concatemersBin[bin2]
                            forZeros = zeroSums[bin2]
                    elif concatemersBin[bin2] == 0:
                        divider = concatemersBin[bin1]
                        forZeros = zeroSums[bin1]
                    else:
                        divider = concatemersBin[bin1] + concatemersBin[bin2]
                        forZeros = zeroSums[bin1] + zeroSums[bin2]

                    #divider = float(concatemersBin[bin1] + concatemersBin[bin2])

                    if divider == 0:
                        divider = 1
                    # if both are zero 
                    norm_data[bin1, bin2] = (hic_data[bin1, bin2] * divider) / float(forZeros)
                    norm_data[bin2, bin1] = norm_data[bin1, bin2]

    elif method == 'concatSum':
        for nbin1, bin1 in enumerate(regRange):
            for bin2 in regRange[nbin1:]:
                if bin1 == bin2:
                    norm_data[bin1, bin2] = 0
                elif hic_data[bin1, bin2] <= mininter:
                    norm_data[bin1, bin2] = 0
                elif bin1 in viewPoint and not bin2 in viewPoint:
                    divider = float(concatemersBin[bin1]*2)
                    norm_data[bin1, bin2] = hic_data[bin1, bin2] / float(divider)
                    norm_data[bin2, bin1] = norm_data[bin1, bin2]
                elif bin2 in viewPoint and not bin1 in viewPoint:
                    divider = float(concatemersBin[bin2]*2)
                    norm_data[bin1, bin2] = hic_data[bin1, bin2] / float(divider)
                    norm_data[bin2, bin1] = norm_data[bin1, bin2]
                else:
                    # get divider
                    if concatemersBin[bin1] == 0:
                        if concatemersBin[bin2] == 0:
                            divider = 1
                        else:
                            divider = concatemersBin[bin2]
                    elif concatemersBin[bin2] == 0:
                        divider = concatemersBin[bin1]
                    else:
                        divider = concatemersBin[bin1] + concatemersBin[bin2]

                    #divider = float(concatemersBin[bin1] + concatemersBin[bin2])

                    if divider == 0:
                        divider = 1
                    # if both are zero 
                    norm_data[bin1, bin2] = hic_data[bin1, bin2] / float(divider)
                    norm_data[bin2, bin1] = norm_data[bin1, bin2]
                    
    # Will remove from divider concatemers counted twice
    elif method == 'concatSum2':
        for nbin1, bin1 in enumerate(regRange):
            for bin2 in regRange[nbin1:]:
                if bin1 == bin2:
                    norm_data[bin1, bin2] = 0
                elif hic_data[bin1, bin2] <= mininter:
                    norm_data[bin1, bin2] = 0
                elif bin1 in viewPoint and not bin2 in viewPoint:
                    divider = float(concatemersBin[bin1]*2)
                    norm_data[bin1, bin2] = hic_data[bin1, bin2] / float(divider)
                    norm_data[bin2, bin1] = norm_data[bin1, bin2]
                elif bin2 in viewPoint and not bin1 in viewPoint:
                    divider = float(concatemersBin[bin2]*2)
                    norm_data[bin1, bin2] = hic_data[bin1, bin2] / float(divider)
                    norm_data[bin2, bin1] = norm_data[bin1, bin2]
                else:
                    # get divider
                    if concatemersBin[bin1] == 0:
                        if concatemersBin[bin2] == 0:
                            divider = 1
                        else:
                            divider = concatemersBin[bin2]
                    elif concatemersBin[bin2] == 0:
                        divider = concatemersBin[bin1]
                    else:
                        divider = concatemersBin[bin1] + concatemersBin[bin2] - hic_data[bin1, bin2]

                    #divider = float(concatemersBin[bin1] + concatemersBin[bin2])

                    if divider == 0:
                        divider = 1
                    # if both are zero 
                    norm_data[bin1, bin2] = hic_data[bin1, bin2] / float(divider)
                    norm_data[bin2, bin1] = norm_data[bin1, bin2]
                    
    elif method == 'concatMtl':
        for nbin1, bin1 in enumerate(regRange):
            for bin2 in regRange[nbin1:]:
                if bin1 == bin2:
                    norm_data[bin1, bin2] = 0
                elif hic_data[bin1, bin2] <= mininter:
                    norm_data[bin1, bin2] = 0
                elif bin1 in viewPoint and not bin2 in viewPoint:
                    divider = float(concatemersBin[bin1]**2)
                    norm_data[bin1, bin2] = hic_data[bin1, bin2] / float(divider)
                    norm_data[bin2, bin1] = norm_data[bin1, bin2]
                elif bin2 in viewPoint and not bin1 in viewPoint:
                    divider = float(concatemersBin[bin2]**2)
                    norm_data[bin1, bin2] = hic_data[bin1, bin2] / float(divider)
                    norm_data[bin2, bin1] = norm_data[bin1, bin2]
                else:
                    # get divider
                    if concatemersBin[bin1] == 0:
                        if concatemersBin[bin2] == 0:
                            divider = 1
                        else:
                            divider = concatemersBin[bin2]
                    elif concatemersBin[bin2] == 0:
                        divider = concatemersBin[bin1]
                    else:
                        divider = concatemersBin[bin1] * concatemersBin[bin2]


                    if divider == 0:
                        divider = 1
                    # if both are zero 
                    norm_data[bin1, bin2] = hic_data[bin1, bin2] / float(divider)
                    norm_data[bin2, bin1] = norm_data[bin1, bin2]
           
    elif method == 'concatSameSum':
        for nbin1, bin1 in enumerate(regRange):
            for bin2 in regRange[nbin1:]:
                if bin1 == bin2:
                    norm_data[bin1, bin2] = 0
                elif hic_data[bin1, bin2] <= mininter:
                    norm_data[bin1, bin2] = 0
                
                else:
                    # get divider
                    if concatemersBin[bin1] == 0:
                        if concatemersBin[bin2] == 0:
                            divider = 1
                        else:
                            divider = concatemersBin[bin2]
                    elif concatemersBin[bin2] == 0:
                        divider = concatemersBin[bin1]
                    else:
                        divider = concatemersBin[bin1] + concatemersBin[bin2]

                    #divider = float(concatemersBin[bin1] + concatemersBin[bin2])

                    if divider == 0:
                        divider = 1
                    # if both are zero 
                    norm_data[bin1, bin2] = hic_data[bin1, bin2] / float(divider)
                    norm_data[bin2, bin1] = norm_data[bin1, bin2]
                    
    # Will remove from divider concatemers counted twice
    elif method == 'concatSameSum2':
        if keep == False:
            for nbin1, bin1 in enumerate(regRange):
                for bin2 in regRange[nbin1:]:
                    # If diagonal or bellow sed interactions
                    if bin1 == bin2 or hic_data[bin1, bin2] <= mininter:
                        pass  # Leave it as zero

                    else:
                        # get divider
                        if concatemersBin[bin1] == 0:
                            if concatemersBin[bin2] == 0:
                                divider = 1
                            else:
                                divider = concatemersBin[bin2]
                        elif concatemersBin[bin2] == 0:
                            divider = concatemersBin[bin1]
                        else:
                            divider = concatemersBin[bin1] + concatemersBin[bin2] - hic_data[bin1, bin2]

                        #divider = float(concatemersBin[bin1] + concatemersBin[bin2])

                        if divider == 0:
                            divider = 1
                        # if both are zero 
                        norm_data[bin1, bin2] = hic_data[bin1, bin2] / float(divider)
                        norm_data[bin2, bin1] = norm_data[bin1, bin2]
                    
        else:
            for ke in keep:
                # If diagonal or bellow sed interactions
                if (ke[0] == ke[1] or hic_data[ke[0], ke[1]] <= mininter or
                    (ke[0] not in regRange or ke[1] not in regRange)):
                    pass  # Leave it as zero
                else:
                    # get divider
                    if concatemersBin[ke[0]] == 0:
                        if concatemersBin[ke[1]] == 0:
                            divider = 1
                        else:
                            divider = concatemersBin[ke[1]]
                    elif concatemersBin[ke[1]] == 0:
                        divider = concatemersBin[ke[0]]
                    else:
                        divider = concatemersBin[ke[0]] + concatemersBin[ke[1]] - hic_data[ke[0], ke[1]]

                    #divider = float(concatemersBin[bin1] + concatemersBin[bin2])

                    if divider == 0:
                        divider = 1
                    # if both are zero 
                    norm_data[ke[0], ke[1]] = hic_data[ke[0], ke[1]] / float(divider)
                    norm_data[ke[1], ke[0]] = norm_data[ke[0], ke[1]]
                    
    elif method == 'justKeep':
        for ke in keep:
            # If diagonal or bellow sed interactions
            if (ke[0] == ke[1] or hic_data[ke[0], ke[1]] <= mininter or
                (ke[0] not in regRange or ke[1] not in regRange)):
                pass  # Leave it as zero
            else:
                norm_data[ke[0], ke[1]] = hic_data[ke[0], ke[1]] 
                norm_data[ke[1], ke[0]] = norm_data[ke[0], ke[1]]
        
    elif method == 'concatSameMlt':
        for nbin1, bin1 in enumerate(regRange):
            for bin2 in regRange[nbin1:]:
                if bin1 == bin2:
                    norm_data[bin1, bin2] = 0
                elif hic_data[bin1, bin2] <= mininter:
                    norm_data[bin1, bin2] = 0
                
                else:
                    # get divider
                    if concatemersBin[bin1] == 0:
                        if concatemersBin[bin2] == 0:
                            divider = 1
                        else:
                            divider = concatemersBin[bin2]
                    elif concatemersBin[bin2] == 0:
                        divider = concatemersBin[bin1]
                    else:
                        divider = concatemersBin[bin1] * concatemersBin[bin2]

                    #divider = float(concatemersBin[bin1] + concatemersBin[bin2])

                    if divider == 0:
                        divider = 1
                    # if both are zero 
                    norm_data[bin1, bin2] = hic_data[bin1, bin2] / float(divider)
                    norm_data[bin2, bin1] = norm_data[bin1, bin2]
                    
    elif method == 'rawNV':
        for nbin1, bin1 in enumerate(regRange):
            for bin2 in regRange[nbin1:]:
                if bin1 == bin2:
                    norm_data[bin1, bin2] = 0
                elif hic_data[bin1, bin2] <= mininter:
                    norm_data[bin1, bin2] = 0
                elif bin1 in viewPoint or bin2 in viewPoint:
                    norm_data[bin1, bin2] = 0
                    norm_data[bin2, bin1] = norm_data[bin1, bin2]
                
                else:
                    norm_data[bin1, bin2] = hic_data[bin1, bin2] 
                    norm_data[bin2, bin1] = norm_data[bin1, bin2]

    # Remove interactions between viewPoints
    if len(viewPoint) > 1:
        comb = list(itertools.combinations(viewPoint, 2))
        for c in comb:
            norm_data[c[0], c[1]] = 0
            norm_data[c[1], c[0]] = 0
            
    return norm_data, concatemersBin

# this is the old version like in TADbit, with a huge list of bins as dict_sect
#def checkOverlapCount(cr1, p1, p2, split_by_bin, resolution, dict_sec):
#    '''
#    Function to get fragment bin location of all the bins included in a fragment
#     and having more than split_by_bin nucleotides inside a bin. Function is called
#     in load_hic_data_from_reads when split_by_bin > 1 and != False
#    '''
#    interacting = []
#    # if first fragments passes thresshold of minimum allowed bp till bin change
#    if (((p1/resolution) + 1) * resolution) - p1 >= split_by_bin:
#        try:
#            interacting.append(dict_sec[(cr1, p1 / resolution)])
#        except:
#            interacting.append(p1 / resolution)
#            
#    # same for second and fragment (we add all posible bins between extremes)
#    if p2 - ((p2/resolution) * resolution) >= split_by_bin:
#        try:
#            interacting += range(dict_sec[(cr1, p1 / resolution)] + 1, dict_sec[(cr1, p2 / resolution)] + 1)
#        except:
#            interacting += range((p1 / resolution) + 1, (p2 / resolution) + 1)
#            
#    # Just in case there is a middle fragment and the last one didnt pass the thresshold
#    else:
#        try:
#            interacting += range(dict_sec[(cr1, p1 / resolution)] + 1, dict_sec[(cr1, p2 / resolution)])
#        except:
#            interacting += range((p1 / resolution) + 1, (p2 / resolution))
#            
#    return interacting

def checkOverlapCount(cr1, p1, p2, split_by_bin, resolution, dict_sec):
    '''
    Function to get fragment bin location of all the bins included in a fragment
     and having more than split_by_bin nucleotides inside a bin. Function is called
     in load_hic_data_from_reads when split_by_bin > 1 and != False
    '''
    interacting = []
    # if first fragments passes thresshold of minimum allowed bp till bin change
    if (((p1/resolution) + 1) * resolution) - p1 >= split_by_bin:
        try:
            interacting.append(dict_sec[cr1][0] + (p1 / resolution))
        except:
            interacting.append(p1 / resolution)

    # same for second and middle fragment (we add all posible bins between extremes)
    if p2 - ((p2/resolution) * resolution) >= split_by_bin:
        try:
            interacting += range(dict_sec[(cr1, p1 / resolution)] + 1, dict_sec[cr1][0] + (p1 / resolution) + 1)
        except:
            interacting += range((p1 / resolution) + 1, (p2 / resolution) + 1)

    # Just in case there is a middle fragment and the last one didnt pass the thresshold
    else:
        try:
            interacting += range(dict_sec[cr1][0] + (p1 / resolution) + 1, dict_sec[cr1][0] + (p1 / resolution))
        except:
            interacting += range((p1 / resolution) + 1, (p2 / resolution))

    return interacting

def checkOverlapPerc(cr1, p1, p2, split_by_bin, resolution, dict_sec):
    '''
    Function to get fragment bin location of all the bins included in a fragment
     and having a proportion greater or equal to split_by_bin inside a bin. 
     Function is called in load_hic_data_from_reads when split_by_bin <= 1 and
     != False
    '''
    interacting = []
    # if first fragments passes thresshold of minimum allowed bp till bin change
    if ((((p1/resolution) + 1) * resolution) - p1) / float(resolution) >= split_by_bin:
        try:
            interacting.append(dict_sec[(cr1, p1 / resolution)])
        except:
            interacting.append(p1 / resolution)
    # same for second fragment
    if (p2 - ((p2/resolution) * resolution)) / float(resolution) >= split_by_bin:
        try:
            interacting += range(dict_sec[(cr1, p1 / resolution)] + 1, dict_sec[(cr1, p2 / resolution)] + 1)
        except:
            interacting += range((p1 / resolution) + 1, (p2 / resolution) + 1)
            
    # Just in case there is a middle fragment and the last one didnt pass the thresshold
    else:
        try:
            interacting += range(dict_sec[(cr1, p1 / resolution)] + 1, dict_sec[(cr1, p2 / resolution)])
        except:
            interacting += range((p1 / resolution) + 1, (p2 / resolution))
    
    return interacting

    
# Look for all the combinations of multiContacts from 0
#to given value
def lookCombiDefinedRange(findMulti, concatTemp, multiGroups):
    # Get all groups of multiContacts
    for sg in range(3, findMulti + 1):
        # create groups of sg multicontacts
        groups = itertools.combinations(sorted(concatTemp), sg)
        # Iterate over each group
        for gr in groups:
            multiGroups[sg][gr] += 1
        
    return multiGroups


# Look for all the combinations of multiContacts for given value
def lookCombiDefined(findMulti, concatTemp, multiGroups):
    # create groups of sg multicontacts
    groups = itertools.combinations(sorted(concatTemp), findMulti)
    # Iterate over each group
    for gr in groups:
        multiGroups[findMulti][gr] += 1
        
    return multiGroups

# Look for all posible combinations of multiContacts in
#a concatemer
def lookCombiAll(findMulti, concatTemp, multiGroups):
    # Get all groups of multiContacts
    for sg in range(3, len(concatTemp) + 1):
        # create groups of sg multicontacts
        groups = itertools.combinations(sorted(concatTemp), sg)
        # Iterate over each group
        for gr in groups:
            multiGroups[sg][gr] += 1
        
    return multiGroups

# Obtain multiContacts from file
def goThroughConcatemerFile(hic_data, line, multiGroups, concatemers,
                                findMulti, resol, nConcat=0,
                                concatTemp=[], prev='',
                                lookComb=lookCombiDefined,
                                frag1Pos=2, lenFrag1Pos=4,
                                frag2Pos=8, lenFrag2Pos=10,
                                chrom1Pos = 1, chrom2Pos=7,
                                RE1Pos=[5,6], RE2Pos=[11,12],
                                split_by_bin=False, seen=[],
                                checkOverlap = checkOverlapCount):
    '''
    :param lookCombiDefined lookComb: Wether you want to retrieve just
        the specified multicontactacts in findMulti (lookCombiDefined),
        all the multicontacts in the range from 3 to findMulti 
        (lookCombiDefinedRange), or all the existing ones (lookCombiAll)
    '''
    line = line.split()
    # If we are in a new concatemer or first one 
    prev_ = line[0].split('#',1)[0]
    if prev != prev_:
        # You can iterate over an empty set if
        #in first fragment
        for k in concatTemp:
            concatemers[k] += 1

        # Look for all combinations of multiContacts
        #if no one specified
        multiGroups = lookComb(findMulti, concatTemp, multiGroups)
        
        # Reset variable
        concatTemp = []
        prev = prev_

        # New concatemer seen
        nConcat += 1
        
        # reset list of seen RFs
        seen = []
        


    # store each apparition of a fragment in a concatemer
    RE11 = line[RE1Pos[0]]
    RE12 = line[RE1Pos[1]]
    
    RE21 = line[RE2Pos[0]]
    RE22 = line[RE2Pos[1]]
    if split_by_bin == False:
        #we store the bin of the mapping start position
        if (RE11, RE12) not in seen:
            concatTemp.append((int(line[frag1Pos]) / resol) + hic_data.section_pos[line[chrom1Pos]][0])
            seen += [(RE11, RE12)]
            
        if (RE21, RE22) not in seen:
            concatTemp.append((int(line[frag2Pos]) / resol) + hic_data.section_pos[line[chrom2Pos]][0])
            seen += [(RE21, RE22)]
        
    
    else:
        # Since fragments from same concatemer are repeated, we could be overcounting the splitting
        #of the same fragment many times, so need a control to count each fragment just once
        # First fragment
        if (RE11, RE12) not in seen:
            ps1 = int(line[frag1Pos])
            ps1_1 = ps1 / resol
            ps12 = ps1 + int(line[lenFrag1Pos])
            ps1_2 = ps12 / resol
            interacting1 = []

            if ps1_1 != ps1_2:
                interacting1 += checkOverlap(line[chrom1Pos], ps1, ps12, split_by_bin, 
                                                 resol, hic_data.sections)
            else:
                try:
                    interacting1 += [hic_data.sections[(line[chrom1Pos], ps1 / resol)]]
                except:
                    interacting1 += [ps1_1]
                    
            seen += [(RE11, RE12)]
            concatTemp += interacting1
                
        # second fragment 
        if (RE21, RE22) not in seen:
            ps2 = int(line[frag2Pos])
            ps2_1 = ps2 / resol
            ps22 = ps2 + int(line[lenFrag2Pos])
            ps2_2 = ps22 / resol
            interacting2 = []

            if ps2_1 != ps2_2:
                interacting2 += checkOverlap(line[chrom2Pos], ps2, ps22, split_by_bin, 
                                                 resol, hic_data.sections)
            else:
                try:
                    interacting2 += [hic_data.sections[(line[chrom2Pos], ps2 / resol)]]
                except:
                    interacting2 += [ps2_1]
                    
            seen += [(RE21, RE22)]
            concatTemp += interacting2

    return multiGroups, concatemers, nConcat, concatTemp, prev, seen

# Obtain multiContacts from file just in a pairwise manner
def goThroughConcatemerFilePairwise(hic_data, line, concatemers,
                                resol, nConcat=0,
                                concatTemp=[], prev='',
                                    concatIdPos=0,
                                frag1Pos=2, lenFrag1Pos=4,
                                frag2Pos=8, lenFrag2Pos=10,
                                chrom1Pos = 1, chrom2Pos=7,
                                RE1Pos=[5,6], RE2Pos=[11,12],
                                split_by_bin=False, seen={},
                                checkOverlap = checkOverlapCount):
    '''
    :param 2 frag1Pos: Index indicating in wich column of the tsv file are
        located the genomic coordinates of the first fragment in the 
        interaction
    :param 8 frag2Pos: Index indicating in wich column of the tsv file are
        located the genomic coordinates of the second fragment in the 
        interaction
    :param False split_by_bin: Set to i) float between 0 to 1 or to ii) integer from 2
        to above to extend to multiple bins the mapped reads cover more than one bin. 
        If i) will include interaction in bin if x% of bin is covered (dangerous since 
        interactions might be lost), if ii) will include interaction in bin if x 
        nucleotides are inside bin.
    '''
    line = line.split()
    # If we are in a new concatemer or first one 
    prev_ = line[0].split('#',1)[0]
    
    #### Este control ya no haria falta creo, para para la prueba dejamos todo igual
    # con la nueva version no haria falta mirar en que concatemer estamos
    if prev != prev_:
        # You can iterate over an empty set if
        #in first fragment
        for k in concatTemp:
            concatemers[k[0]].add(k[1])

        # Reset variable
        concatTemp = []
        prev = prev_

        # New concatemer seen
        nConcat += 1
        
        # reset list of seen RFs
        seen = {}


    # store each apparition of a fragment in a concatemer
    RE11 = line[RE1Pos[0]]
    RE12 = line[RE1Pos[1]]
    
    RE21 = line[RE2Pos[0]]
    RE22 = line[RE2Pos[1]]
    # store each apparition of a fragment in a concatemer
    concatID = line[concatIdPos].split('#')[0]
    if split_by_bin == False:
        #we store the bin of the mapping start position
        if (RE11, RE12) not in seen:
            pos = (int(line[frag1Pos]) / resol) + hic_data.section_pos[line[chrom1Pos]][0]
            concatTemp.append([pos, concatID])
            seen[(RE11, RE12)] = pos
        else:
            concatTemp.append([seen[(RE11, RE12)], concatID])
            
        if (RE21, RE22) not in seen:
            pos = (int(line[frag2Pos]) / resol) + hic_data.section_pos[line[chrom2Pos]][0]
            concatTemp.append([pos,concatID])
            seen[(RE21, RE22)] = pos
        else:
            concatTemp.append([seen[(RE21, RE22)], concatID])
    else:
        # Since fragments from same concatemer are repeated, we could be overcounting the splitting
        #of the same fragment many times, so need a control to count each fragment just once
        # First fragment
        if (RE11, RE12) not in seen:
            ps1 = int(line[frag1Pos])
            ps1_1 = ps1 / resol
            ps12 = ps1 + int(line[lenFrag1Pos])
            ps1_2 = ps12 / resol
            interacting1 = []

            if ps1_1 != ps1_2:
                interacting1 += checkOverlap(line[chrom1Pos], ps1, ps12, split_by_bin, 
                                                 resol, hic_data.sections)
            else:
                try:
                    interacting1 += [hic_data.sections[(line[chrom1Pos], ps1 / resol)]]
                except:
                    interacting1 += [ps1_1]
                    
            seen[(RE11, RE12)] = interacting1
            concatTemp += [[interacting1, concatID]]
        else:
            concatTemp += [[seen[(RE11, RE12)], concatID]]
            
        # second fragment 
        if (RE21, RE22) not in seen:
            ps2 = int(line[frag2Pos])
            ps2_1 = ps2 / resol
            ps22 = ps2 + int(line[lenFrag2Pos])
            ps2_2 = ps22 / resol
            interacting2 = []

            if ps2_1 != ps2_2:
                interacting2 += checkOverlap(line[chrom2Pos], ps2, ps22, split_by_bin, 
                                                 resol, hic_data.sections)
            else:
                try:
                    interacting2 += [hic_data.sections[(line[chrom2Pos], ps2 / resol)]]
                except:
                    interacting2 += [ps2_1]
                    
            seen[(RE21, RE22)] = interacting2
            concatTemp += [[interacting2, concatID]]
            
        else:
            concatTemp += [[seen[(RE21, RE22)], concatID]]

    return concatemers, nConcat, concatTemp, prev, seen
    

# Open tsv and obtain multiContact frecuencies
def getMultiAndConcatemersPerBin(hic_data, tsvFile, resol, locusCh=False,
                                 regRange=False, returnNconcat = False,
                                findMulti=False, lookComb=lookCombiDefined,
                                split_by_bin=False):
    '''
    Function to get the number of concatemers were a bin of interest is 
        appearing and the multi contact groups (is bin based, so all 
        fragments which start inside of a bin margin will be joined). 
        Assumes integer based chromosomes, where mt would be 26 in humans
    :param False returnNconcat: wether you want or not nConcat to be
        returned
    :param False findMulti: Integer if you want to look just for multi
        contacts with that amount o members
    :param lookCombiDefined lookComb: Wether you want to retrieve just
        the specified multicontactacts in findMulti (lookCombiDefined),
        all the multicontacts in the range from 3 to findMulti 
        (lookCombiDefinedRange), or all the existing ones (lookCombiAll)
    :param False split_by_bin: Set to i) float between 0 to 1 or to ii) integer from 2
        to above to extend to multiple bins the mapped reads cover more than one bin. 
        If i) will include interaction in bin if x% of bin is covered (dangerous since 
        interactions might be lost), if ii) will include interaction in bin if x 
        nucleotides are inside bin.
    '''
    # Prepare multi contact retrieving stratey
    if findMulti == False:
        lookComb=lookCombiAll
        
        # Variable to store multiContact groups
        multiGroups = {}
        # asume maximum of 20 multiContacts to avoid if statments
        for sg in range(3, 21):
            multiGroups[sg] = defaultdict(int)
        
    else:
        multiGroups = {findMulti:defaultdict(int)}    
 
    # variable to store id
    prev = ''
    concatemers = defaultdict(int)

    # Use set to remove duplicates or fragments from same bin in a concatemer
    concatTemp = []
    nConcat = 0
    
    if split_by_bin <= 1:
        checkOverlap = checkOverlapPerc
    else:
        checkOverlap = checkOverlapCount

    

    with open(tsvFile, 'r') as f:
        # Skip initial comments that starts with #
        while True:
            line = f.readline()
            # break while statement if it is not a comment line
            # i.e. does not startwith #
            if not line.startswith('#'):
                break
        # Run current line (first one)
        multiGroups, concatemers, nConcat, concatTemp, prev, seen = goThroughConcatemerFile(hic_data, 
                                line, multiGroups, concatemers, findMulti, resol, nConcat=0,
                                concatTemp=[], prev='', lookComb=lookComb,split_by_bin=split_by_bin, 
                                                checkOverlap = checkOverlap)
        # Go for next lines
        for line in f:
            multiGroups, concatemers, nConcat, concatTemp, prev, seen = goThroughConcatemerFile(hic_data, line, multiGroups, concatemers,
                                findMulti, resol, nConcat=nConcat,
                                concatTemp=concatTemp, prev=prev, lookComb=lookComb,split_by_bin=split_by_bin, 
                                                checkOverlap = checkOverlap, seen=seen)

    # Add last concatemer of file
    for k in concatTemp:
        concatemers[k] += 1
    multiGroups = lookComb(findMulti, concatTemp, multiGroups)

    
        

    # Get genomic coordinates for region of interest
    if locusCh == True:
        regionStart = min(regRange)
        regionEnd = max(regRange) 

        ## modify if we want data from all genome
        # Get all the fragments that start inside this coordinates
        keys = [k for k in concatemers.keys() if (regionStart <= 
                                                            k <= 
                                                            regionEnd)]
        regConcatemers = defaultdict(int)
        for k in keys:
            regConcatemers[k] += concatemers[k]
            
    # Or not
    else:
        regConcatemers = concatemers

                

    ##
    # Remove multiGroups keys with no data
    keys = multiGroups.keys()
    for sg in keys:
        if len(multiGroups[sg]) == 0:
            del multiGroups[sg]

    if returnNconcat == False:
        return regConcatemers, multiGroups
    else:
        return regConcatemers, multiGroups, nConcat
    
    

    

# Open tsv and obtain multiContact frecuencies in a pairwise manner
def getConcatemersPerBin(hic_data, tsvFile, resol, locusCh=False,
                         regRange = False, returnNconcat = False,
                         split_by_bin=False):
    '''
    Function to get the number of concatemers were a bin of interest is 
        appearing (is bin based, so all fragments which start inside of
        a bin margin will be joined)
    :param False split_by_bin: Set to i) float between 0 to 1 or to ii) integer from 2
        to above to extend to multiple bins the mapped reads cover more than one bin. 
        If i) will include interaction in bin if x% of bin is covered (dangerous since 
        interactions might be lost), if ii) will include interaction in bin if x 
        nucleotides are inside bin.
        
    :param False returnNconcat: wether you want or not nConcat to be
        returned
    '''
    # variable to store id
    prev = ''
    concatemers = defaultdict(int)
    
    concatTemp = []
    nConcat = 0
    
    if split_by_bin <= 1:
        checkOverlap = checkOverlapPerc
    else:
        checkOverlap = checkOverlapCount

    with open(tsvFile, 'r') as f:
        # Skip initial comments that starts with #
        while True:
            line = f.readline()
            if not line.startswith('#'):
                break
        # Run current line (first one)
        concatemers, nConcat, concatTemp, prev, seen = goThroughConcatemerFilePairwise(hic_data, 
                                                line, concatemers, resol, nConcat=0, 
                                                concatTemp=[], prev='', split_by_bin=split_by_bin, 
                                                checkOverlap = checkOverlap)
        
        # Go for next lines
        for line in f:
            concatemers, nConcat, concatTemp, prev, seen = goThroughConcatemerFilePairwise(hic_data, 
                                                line, concatemers, resol, nConcat=nConcat,
                                                concatTemp=concatTemp, prev=prev, 
                                                split_by_bin=split_by_bin, checkOverlap = checkOverlap,
                                                seen=seen)

    # Add last concatemer of file
    for k in concatTemp:
        concatemers[k] += 1


     # filter genomic coordinates outside region of interest
    if locusCh == True:
        regionStart = min(regRange)
        regionEnd = max(regRange) 

        ## modify if we want data from all genome
        # Get all the fragments that start inside this coordinates
        keys = [k for k in concatemers.keys() if (regionStart <= 
                                                            k <= 
                                                            regionEnd)]
        regConcatemers = defaultdict(int)
        for k in keys:
            regConcatemers[k] += concatemers[k]
            
    # Or not
    else:
        regConcatemers = concatemers
            

    if returnNconcat == False:
        return regConcatemers
    else:
        return regConcatemers, nConcat

# Prepare randomised shufling of multicontacts
def randomiseMultiGroups(multiGroups, multiLevel):
    # we get all keys from the multiGroups
    groups = []
    for k in multiGroups[multiLevel].keys():
        for i in range(multiGroups[multiLevel][k]):
            groups.append(k)

    # And split each interacting bin
    if multiLevel != len(groups[0]):
        print 'WARNING: multiLevel different from obtained groups size'
    allKeys = []
    for k in groups:
        for i in range(multiLevel):
            allKeys.append(k[i])

    # Then shuffle them
    random.shuffle(allKeys)

    # Now merge them in groups of multiLevel size
    groups = {multiLevel: defaultdict(int)}
    for i in range(0, len(allKeys), 3):
        groups[multiLevel][tuple(allKeys[i:i+multiLevel])] += 1

    return groups

# Function to normalise by frecuencies
def MultiNorm(hic_data, regionStartBin, regionEndBin, tsvFile, resol, locusCh=False,
              method='', multiLevel=2, zeroPercent=False, mininter=0, keep=False, 
              returnNconcat=False, random=False, multResult=100, split_by_bin=False):
    '''
    :param regionStartBin: integer indicating the bin FROM wich we will store 
        normalised data in the new hic_data object
    :param regionEndBin: integer indicating the bin TO wich (itself included)
        we will store normalised data in the new hic_data object
    :param False keep: if you have a list of bin coordinates that should be the only 
        ones to be kept. e.j. [(1,3), (235,76), ...]
    :param 2 multiLevel: To normalise taking into account n-wise interactions
    :param False locusCh: Add chromosome name if focus in a region of interest
    :param 100 multResult: integer or floar with wich multiply the resulting
        interaction scores. Default as 100 since the result is like a percentaje
    :param False random: Set to true if you want to normalise random distributions
        of the multicontacts
    :param False split_by_bin: Set to i) float between 0 to 1 or to ii) integer from 2
        to above to extend to multiple bins the mapped reads cover more than one bin. 
        If i) will include interaction in bin if x% of bin is covered (dangerous since 
        interactions might be lost), if ii) will include interaction in bin if x 
        nucleotides are inside bin.
    '''
    
    sumStart = regionStartBin
    sumEnd = regionEndBin
    # Select till which point are we going to take into account interactions at
    #the time to normalise
    # Get iteration bins ranges
    regRange = set(range(regionStartBin, regionEndBin + 1))
    normRange = set(range(sumStart, sumEnd + 1))
    
    # Get concatemer appearances per bin
    if multiLevel == 2:
        if returnNconcat == False:
            concatemersBin = getConcatemersPerBin(hic_data, tsvFile, resol, 
                                                  locusCh, regRange, returnNconcat,
                                                  split_by_bin=split_by_bin)
        else:
            concatemersBin, nConcat = getConcatemersPerBin(hic_data, tsvFile,
                                             resol, locusCh, regRange, returnNconcat,
                                             split_by_bin=split_by_bin)
    # If we are going to check for multi contacts
    elif multiLevel > 2:
        if returnNconcat == False:
            concatemersBin, multiGroups = getMultiAndConcatemersPerBin(hic_data, tsvFile, 
                                                                       resol, locusCh, regRange,
                                                                       returnNconcat,
                                                                       split_by_bin=split_by_bin,
                                                                       findMulti=multiLevel)
        else:
            concatemersBin, multiGroups, nConcat = getMultiAndConcatemersPerBin(hic_data, 
                                                                                tsvFile, resol,
                                                                                locusCh, 
                                                                                regRange,
                                                                                returnNconcat,
                                                                                split_by_bin=split_by_bin,
                                                                                findMulti=multiLevel)
        # Clue to know how many multi contact groups we have
        #print 'Longest multicontact: %s' %max(multiGroups.keys())
    
    # If we want random distributions of multigroups
    if random == True:
        # This procedure has not been optimised for pairwise interactions
        if multiLevel <= 2:
            print 'WARNING: randomisation procedure has not been optimised \
for pairwise interactions. Wont randomise'
        else:
            multiGroups = randomiseMultiGroups(multiGroups, multiLevel)
    
    
    zeroSums = {}
    ## Filter data
    if zeroPercent != False:
        # For that, we turn it into zero
        # If badsThres is given in percentaje
        if zeroPercent <= 1:
            nbins = regionEndBin - regionStartBin + 1
            # Minimum number of non zero values
            badsThres = nbins * zeroPercent
        # If is given in minimum number of interactions (for one use 1.1)
        else:
            badsThres = zeroPercent
        bads = set()

        # First get 
        for bin1 in regRange:
                # If we are in the viewpoint just avoid diagonal
                zeroSums[bin1] = sum(i != bin1 and 
                                     not hic_data[bin1, i]
                                     for i in normRange)
                if zeroSums[bin1] > badsThres:
                    bads.add(bin1)

                # get non zero cells
                suma = sum(i != bin1 and 
                            hic_data[bin1, i]
                             for i in range(sumStart, sumEnd + 1))
                # If we have all zeros leave it as one for the enumerater, since
                #interaction value will alway be zero
                if suma == 0:
                    zeroSums[bin1] = 1

                else:
                    zeroSums[bin1] = ((zeroSums[bin1] + 1)/ float(suma))

        # Remove filtered columns
        regRange = regRange - bads
        normRange = normRange - bads
    
    # Now turn set to list again to recover order
    regRange = list(sorted(regRange))
    normRange = list(sorted(normRange))
     
    if zeroPercent == False:
        # get number of zeros if no filtering step
        for bin1 in regRange:
            # If we are in the viewpoint just avoid diagonal
            zeroSums[bin1] = sum(i != bin1 and 
                                 not hic_data[bin1, i]
                                 for i in normRange)
            # If all points have data no adjustment required, so value = 1
            # pass

            suma = sum(i != bin1 and 
                            hic_data[bin1, i]
                             for i in range(sumStart, sumEnd + 1))
            # If we have all zeros leave it as one for the enumerater, since
            #interaction value will alway be zero
            if suma == 0:
                zeroSums[bin1] = 1

            else:
                zeroSums[bin1] = ((zeroSums[bin1] + 1)/ float(suma))
                #zeroSums[bin1] = zeroSums[bin1]
                
              
    if multiLevel == 2:
        # create HiC data for normalised interactions
        dict_sec = hic_data.sections
        genome_seq = hic_data.chromosomes
        size = sum(genome_seq[crm] for crm in genome_seq)
        norm_data = HiC_data((), size, genome_seq, dict_sec, resolution=resol)
    
        # Will remove from divider concatemers counted twice
        if method == 'concatSameSum2':
            if keep == False:
                for nbin1, bin1 in enumerate(regRange):
                    for bin2 in regRange[nbin1:]:
                        # If diagonal or bellow sed interactions
                        if bin1 == bin2 or hic_data[bin1, bin2] <= mininter:
                            pass  # Leave it as zero

                        else:
                            # get divider
                            #if concatemersBin[bin1] == 0:
                            #    if concatemersBin[bin2] == 0:
                            #        divider = 1
                            #    else:
                            #        divider = concatemersBin[bin2]
                            #elif concatemersBin[bin2] == 0:
                            #    divider = concatemersBin[bin1]
                            #else:
                            divider = concatemersBin[bin1] + concatemersBin[bin2] - hic_data[bin1, bin2]

                            #divider = float(concatemersBin[bin1] + concatemersBin[bin2])

                            #if divider == 0:
                            #    divider = 1
                            # if both are zero 
                            norm_data[bin1, bin2] = (hic_data[bin1, bin2] / float(divider)) * multResult
                            norm_data[bin2, bin1] = norm_data[bin1, bin2]

            else:
                for ke in keep:
                    # If diagonal or bellow sed interactions
                    if (ke[0] == ke[1] or hic_data[ke[0], ke[1]] <= mininter or
                        (ke[0] not in regRange or ke[1] not in regRange)):
                        pass  # Leave it as zero
                    else:
                        if not hic_data[ke[0], ke[1]] <= mininter:
                        # get divider
                        #if concatemersBin[ke[0]] == 0:
                        #    if concatemersBin[ke[1]] == 0:
                        #        divider = 1
                        #    else:
                        #        divider = concatemersBin[ke[1]]
                        #elif concatemersBin[ke[1]] == 0:
                        #    divider = concatemersBin[ke[0]]
                        #else:
                            divider = concatemersBin[ke[0]] + concatemersBin[ke[1]] - hic_data[ke[0], ke[1]]

                        #divider = float(concatemersBin[bin1] + concatemersBin[bin2])

                        #if divider == 0:
                        #    divider = 1
                        # if both are zero 
                        norm_data[ke[0], ke[1]] = hic_data[ke[0], ke[1]] / float(divider)
                        norm_data[ke[1], ke[0]] = norm_data[ke[0], ke[1]]
                        
    elif multiLevel >= 3:
        # Will remove from divider concatemers counted twice
        if method == 'concatSameSum2':
            # Check for existance of that level
            if multiLevel in multiGroups.keys():
                # Get multiContact positions keys
                keys = multiGroups[multiLevel].keys()
                # Select the ones inside our region
                keys = [k for k in keys if (k[0] >= regionStartBin and k[-1] <= regionEndBin)]

                # Create new dict for normalised values
                focusMultiGroups = {}
                # Normalise
                for k in keys:
                    # first get presence of all fragments in the whole dataset
                    #meaning, number of times this fragment appeared
                    #divider = sum(concatemersBin[k[nm]] for nm in range(multiLevel))
                    # then remove concatemers counted more than once, maybe because
                    #two of the fragments from k were present in them
                    # look for appearances of more than one member of the multiContact            
                    #divider -= sum(multiGroups[multiLevel][k2] * (
                                        # this line will give as how many times we counted extra this
                                        #concatemer (given presence of fragments from k)
                    #                    sum(1 for kk in k if kk in k2) - 1)
                    #               for k2 in keys 
                                   # To avoid multiplication by negative numbers, we just let pass the check
                                   #concatemers were we know that at least 2 fragments from k are present
                    #               if sum(1 for kk in k if kk in k2) >= 2)
                    divider = sum(multiGroups[multiLevel][k2]
                               for k2 in keys 
                               if sum(1 for kk in k if kk in k2) >= 1)
                    if divider < 0:
                        print divider, [concatemersBin[k[nm]] for nm in range(multiLevel)], multiGroups[multiLevel][k]

                    focusMultiGroups[k] = (multiGroups[multiLevel][k] / float(divider)) * multResult
                norm_data = focusMultiGroups
                
            # If we dont have that level of multiContacts
            else:
                print 'WWARNING: No multicontacts with %s level' %multiLevel
                norm_data = {}
    
    if returnNconcat == False:      
        return norm_data, concatemersBin
    else:
        return norm_data, concatemersBin, nConcat




# Outborder areas are counted as zeros
def median_filter_dict(data, filter_size, start=False, end=False):
    '''
    Function to apply median filter to a matrix of any number of
        dimension. Is prepared for genome interaction data, where 
        interaction between A and B and B and A are the same, so
        just interactions where A<B are mantained and processed.
        Areas outside the matrix or in A>B are treated as zeros.
        Matrix is introduced as a dictionary
    
    :param data0: dictionary with data to be median filtered
    :param filter_size: n x n dimensiones of the submatrix size to
        apply median filtering. MUST be odd
    '''
    if start == False:
        start = np.amin(data.keys())
    if end == False:
        end = np.amax(data.keys())
        
    # Be sure that filter_size is odd
    if filter_size % 2 == 0:
        print 'WARNING: set filter_size as odd number'
    # Obtain number of dimensions in our dataset
    keylen = len(data.keys()[0])
    
    # Create variable to store output
    data2 = {}
    # Integer with number of bins up and down from each filtered submatrix
    indexer = filter_size // 2
    # Create list of coordinates for the position coordinates (away
    #from filtered bin) to be check given filter_size
    focusRng = range(-indexer, filter_size-indexer)
    window = list(itertools.product(focusRng,repeat=keylen))
    
    # Half bins to be check (integer)
    index = len(window) // 2
    for k in data:
        data2[k] = sorted(
            ## If out of borders in y or x axis
            0 if (
                # If any coordinate before begining point
                min(k[nkl]+a[nkl] for nkl in range(keylen)) < start
                # If any coordinate duplicated
                #we will always use the sorted order fo each N-wise combination
                or [k[nkl]+a[nkl] for nkl in range(keylen)] != sorted(k[nkl]+a[nkl] for nkl in range(keylen))
                # At end
                or sum(end < (k[nkl]+a[nkl]) for nkl in range(keylen)) > 0
            ## If inside borders (beware of empty cells)
            # max between 0 and None to obtain 0
            ) else max(data.get(tuple(k[nkl]+a[nkl] for nkl in range(keylen))), 0)
            for a in window
        # Take median value
        )[index]
    return data2
    

def reNorm(matrix, minFull = 0):
    '''
    Funtion to normalise pairwarise contact matrices in order to remove
        the view-point bias.
    :param matrix: list o list with the 2D interaction data
    :param 0 minFull: Minimum number o cells that must have data in a
        matrix row in order to take it into account (otherwise is turned
        to all zeros)
    '''
    matrix2 = [[0 for i in range(len(matrix))] for ii in range(len(matrix))]
    for i in range(len(matrix)):
        isum = sum(matrix[i])
        for j in range(len(matrix)):
            jsum = sum(matrix[j])
            
            divi = isum + jsum
            if divi == 0:
                divi = 1

            matrix2[i][j] = matrix[i][j] / float(divi - matrix[i][j])

    # Remove ares with low data
    for i in range(len(matrix)):
        full = sum(1 for j in range(len(matrix)) if matrix2[i][j] != 0)
        if full < minFull:
            for j in range(len(matrix)):
                matrix2[i][j] = 0
    return matrix2



def multiReNorm(data):
    '''
    Function to normalise for viewPoint the multi contact data
    :param data: dictionary with multiContact values
    '''

    # Create variable to store output
    data2 = {}

    # get appearances of each position
    interSums = defaultdict(int)
    for da in data:
        for d in da:
            interSums[d] += data[da]

    # Normalise
    for da in data.keys():
        # sum all appearances of fragments in the multicontact group
        divisor = sum(interSums[da[nd]] for nd in range(len(da)))

        # look for appearances of more than one member of the multiContact            
        divisor -= sum(data[da2] * (sum(1 for d in da if d in da2) - 1) 
                       for da2 in data 
                       # These concatemers with more than one member were counted more than once
                        # sum(1 for d in da if d in da2) - 1 indicates how many extra times were counted 
                       if sum(1 for d in da if d in da2) >= 2)

        data2[da] = data[da] / float(divisor)


    return data2
