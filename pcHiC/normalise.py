from collections import defaultdict
from pytadbit                     import HiC_data


# Obtain multiContacts from file just in a pairwise manner
def goThroughConcatemerFilePairwise_noMultiTSV(hic_data, line, concatemers,
                                resol, nConcat=0):
    '''
    Function to obtain interaction frecuencies from normal TSV with no 
        multiContact data spected
    
    '''
    line = line.split()
   
    # store each apparition of a fragment in a concatemer
    #we store the bin of the mapping start position
    fragment1 = (int(line[2]) / resol) + hic_data.section_pos[line[1]][0]
    concatemers[fragment1] += 1
    
    fragment2 = (int(line[8]) / resol) + hic_data.section_pos[line[7]][0]
    concatemers[fragment2] += 1
    
    # New concatemer seen
    nConcat += 1
    
    return concatemers, nConcat

    


    

# Open tsv and obtain multiContact frecuencies in a pairwise manner
def getConcatemersPerBin_noMultiTSV(hic_data, tsvFile, resol, locusCh=False,
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
        concatemers, nConcat = goThroughConcatemerFilePairwise_noMultiTSV(hic_data, line, concatemers,
                                resol, nConcat=0)
        
        # Go for next lines
        for line in f:
            concatemers, nConcat = goThroughConcatemerFilePairwise_noMultiTSV(hic_data, line, concatemers,
                                resol, nConcat=nConcat)

    
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




                                            
# Normalise by frecuencies given the presence of each interacting fragment                                           
def frecuenciesNorm(hic_data, resol, regRange, concatemersBin, multResult=100, 
                    keep=False, mininter=0):
    # create HiC data for normalised interactions
    dict_sec = hic_data.sections
    genome_seq = hic_data.chromosomes
    size = sum(genome_seq[crm] for crm in genome_seq)
    norm_data = HiC_data((), size, genome_seq, dict_sec, resolution=resol)

    # Will remove from divider concatemers counted twice
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
        
    return norm_data
                                         

# Whole Function to normalise by frecuencies
def pcHiCnorm_multiway(hic_data, regionStartBin, regionEndBin, tsvFile, resol, locusCh=False,
              zeroPercent=False, mininter=0, keep=False, 
              returnNconcat=False, multResult=100):
    '''
    :param False keep: if you have a list of bin coordinates that should be the only 
        ones to be kept. e.j. [(1,3), (235,76), ...]
    :param False locusCh: Add chromosome name if focus in a region of interest
    :param 100 multResult: integer or floar with wich multiply the resulting
        interaction scores. Default as 100 since the result is like a percentaje
    '''
    
    sumStart = regionStartBin
    sumEnd = regionEndBin
    # Select till which point are we going to take into account interactions at
    #the time to normalise
    # Get iteration bins ranges
    regRange = set(range(regionStartBin, regionEndBin + 1))
    normRange = set(range(sumStart, sumEnd + 1))
    
    # Get concatemer appearances per bin
    if returnNconcat == False:
        concatemersBin = getConcatemersPerBin(hic_data, tsvFile, resol, 
                                              locusCh, regRange, returnNconcat)
    else:
        concatemersBin, nConcat = getConcatemersPerBin(hic_data, tsvFile,
                                         resol, locusCh, regRange, returnNconcat)


    
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
                
          
    ## Lets normalise
    norm_data = frecuenciesNorm(hic_data, resol, regRange, concatemersBin, multResult=multResult, 
                    keep=keep, mininter=mininter)
    
    
    
    if returnNconcat == False:      
        return norm_data, concatemersBin
    else:
        return norm_data, concatemersBin, nConcat


    
    
