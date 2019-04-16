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
                                concatTemp=set(), prev='',
                           lookComb=lookCombiDefined):
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
        concatTemp = set()
        prev = prev_

        # New concatemer seen
        nConcat += 1


    # store each apparition of a fragment in a concatemer
    #we store the bin of the mapping start position
    # Since its a set, will store aparitions just once
    concatTemp.add((int(line[2]) / resol) + hic_data.section_pos[line[1]][0])
    concatTemp.add((int(line[8]) / resol) + hic_data.section_pos[line[7]][0])

    return multiGroups, concatemers, nConcat, concatTemp, prev

# Obtain multiContacts from file just in a pairwise manner
def goThroughConcatemerFilePairwise(hic_data, line, concatemers,
                                resol, nConcat=0,
                                concatTemp=set(), prev='',
                                frag1Pos=2, frag2Pos=8,
                                chrom1Pos = 1, chrom2Pos=7):
    '''
    :param lookCombiDefined lookComb: Wether you want to retrieve just
        the specified multicontactacts in findMulti (lookCombiDefined),
        all the multicontacts in the range from 3 to findMulti 
        (lookCombiDefinedRange), or all the existing ones (lookCombiAll)
    :param 2 frag1Pos: Index indicating in wich column of the tsv file are
        located the genomic coordinates of the first fragment in the 
        interaction
    :param 8 frag2Pos: Index indicating in wich column of the tsv file are
        located the genomic coordinates of the second fragment in the 
        interaction
    '''
    line = line.split()
    # If we are in a new concatemer or first one 
    prev_ = line[0].split('#',1)[0]
    if prev != prev_:
        # You can iterate over an empty set if
        #in first fragment
        for k in concatTemp:
            concatemers[k] += 1

        # Reset variable
        concatTemp = set()
        prev = prev_

        # New concatemer seen
        nConcat += 1


    # store each apparition of a fragment in a concatemer
    #we store the bin of the mapping start position
    # Since its a set, will store aparitions just once
    concatTemp.add((int(line[frag1Pos]) / resol) + hic_data.section_pos[line[chrom1Pos]][0])
    concatTemp.add((int(line[frag2Pos]) / resol) + hic_data.section_pos[line[chrom2Pos]][0])

    return concatemers, nConcat, concatTemp, prev

    

# Open tsv and obtain multiContact frecuencies
def getMultiAndConcatemersPerBin(hic_data, tsvFile, resol, locusCh=False,
                                 regRange=False, returnNconcat = False,
                                findMulti=False, lookComb=lookCombiDefined):
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
    concatTemp = set()
    nConcat = 0

    

    with open(tsvFile, 'r') as f:
        # Skip initial comments that starts with #
        while True:
            line = f.readline()
            # break while statement if it is not a comment line
            # i.e. does not startwith #
            if not line.startswith('#'):
                break
        # Run current line (first one)
        multiGroups, concatemers, nConcat, concatTemp, prev = goThroughConcatemerFile(hic_data, line, multiGroups, concatemers,
                                findMulti, resol, nConcat=0,
                                concatTemp=set(), prev='', lookComb=lookComb)
        # Go for next lines
        for line in f:
            multiGroups, concatemers, nConcat, concatTemp, prev = goThroughConcatemerFile(hic_data, line, multiGroups, concatemers,
                                findMulti, resol, nConcat=nConcat,
                                concatTemp=concatTemp, prev=prev, lookComb=lookComb)

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
            regConcatemers[k] += concatemers[locusCh][k]
            
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
                         regRange = False, returnNconcat = False):
    '''
    Function to get the number of concatemers were a bin of interest is 
        appearing (is bin based, so all fragments which start inside of
        a bin margin will be joined)
        
    :param False returnNconcat: wether you want or not nConcat to be
        returned
    '''
    # variable to store id
    prev = ''
    concatemers = defaultdict(int)
    
    concatTemp = set()
    nConcat = 0

    with open(tsvFile, 'r') as f:
        # Skip initial comments that starts with #
        while True:
            line = f.readline()
            # break while statement if it is not a comment line
            # i.e. does not startwith #
            if not line.startswith('#'):
                break
        # Run current line (first one)
        concatemers, nConcat, concatTemp, prev = goThroughConcatemerFilePairwise(hic_data, line, concatemers,
                                resol, nConcat=0,
                                concatTemp=set(), prev='')
        
        # Go for next lines
        for line in f:
            concatemers, nConcat, concatTemp, prev = goThroughConcatemerFilePairwise(hic_data, line, concatemers,
                                resol, nConcat=nConcat,
                                concatTemp=concatTemp, prev=prev)

    # Add last concatemer of file
    for k in concatTemp:
        concatemers[k] += 1


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
            regConcatemers[k] += concatemers[locusCh][k]
            
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
              returnNconcat=False, random=False, multResult=100):
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
                                                  locusCh, regRange, returnNconcat)
        else:
            concatemersBin, nConcat = getConcatemersPerBin(hic_data, tsvFile,
                                             resol, locusCh, regRange, returnNconcat)
    # If we are going to check for multi contacts
    elif multiLevel > 2:
        if returnNconcat == False:
            concatemersBin, multiGroups = getMultiAndConcatemersPerBin(hic_data, tsvFile, 
                                                                       resol, locusCh, regRange,
                                                                       returnNconcat)
        else:
            concatemersBin, multiGroups, nConcat = getMultiAndConcatemersPerBin(hic_data, 
                                                                                tsvFile, resol,
                                                                                locusCh, 
                                                                                regRange,
                                                                                returnNconcat)
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
                            norm_data[bin1, bin2] = (hic_data[bin1, bin2] / float(divider)) * multResult
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
                    # first get presence of all fragments 
                    divider = sum(concatemersBin[k[nm]] for nm in range(multiLevel))
                    # then remove duplicated ones
                    # look for appearances of more than one member of the multiContact            
                    divider -= sum(multiGroups[multiLevel][k2] * (sum(1 for kk in k if kk in k2) - 1) 
                                   for k2 in keys 
                                   # These concatemers with more than one member were counted more than once
                                    # sum(1 for d in da if d in da2) - 1 indicates how many extra times were counted 
                                   if sum(1 for kk in k if kk in k2) >= 2)
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


    
    

