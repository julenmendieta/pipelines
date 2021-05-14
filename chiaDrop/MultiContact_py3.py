import collections
from collections import defaultdict
from pytadbit                     import HiC_data
import itertools
import copy
import os
import h5py
import numpy as np
import random
from pytadbit.parsers.hic_bam_parser import filters_to_bin
from pysam                           import AlignmentFile
from collections                     import OrderedDict



#  Function to write the TSV files from chiadrop files
def fromCHIADROPtoTSV_chiaPaper(files, outPath, all_chromLengths,
                    filterQual=21,
                    hic_data_style=True, chrCode = '',
                    removeSingle=True):

    '''
    Function to generate TSV files from the ChIA-drop ones
    :param files: list with paths in which the input ChIA-drop files are located
    :param outPath: String with path in which to store the output tsv files
    :param all_chromLengths: Dictionary with reference genomeas key and
        a text including the lengts of each chromosome like in SAM format
    :param 21 filterQual: Minimum quality accepted
    :param '' chrCode: can be chr for example, in cases in which process data does
        not include it
    :param True hic_data_style: Wether to write the tsv in hic_data style readble
        by TADbit, or just keep minimum posible interactions to be used by our
        pipeline. The first one basically shows all vs all interactions for each
        concatemer, whereas the second one shows the first fragment interaction
        against all the others, saving some space
    :param True removeSingle: Set to True in order to remove singleton fragments
    '''

    # Set it to de text you want to add before eery chromomse.
    nConcatDict = {}
    #doubleRE = {}
    for filepath in sorted(files):
        cmd = ''
        if filepath.endswith('txt'):
            id1 = ''.join(filepath.split('/')[-1].split('.')[0].split('ChIA-Drop-'))
            print('## ' + id1 + ' ##')

            if hic_data_style:
                outfile = outPath + '%s_full.tsv' %(id1)
            else:
                outfile = outPath + '%s_short.tsv' %(id1)

            #viewPoint = proposedView[id1]

            # get reference genome chromosome coordinates for this sample
            chromLengths = all_chromLengths

            # Create list with concatemers with repeated RE
            #doubleRE[id1] = {}

            # Keep record of number of concatemers check
            nConcat = 0
            nConcatDict[id1] = {}
            # Keep record of number of concatermers with duplicated RE
            #nConcatDup = 0

            ##  First to write, the header of the tsv
            for h in chromLengths.split('\n'):
                if len(h) != 0:
                    cmd += '# CRM %s\n' %'\t'.join(h.split())
            #  Create file and write first lines
            with open(outfile,"w") as f:
                f.write(cmd)

            ##  Then write the data
            # start adding the contacts
            concatemer = []

            # index of positions of interest
            FragsPos = 4

            # strand position always in forward
            # I believe this should mean forward?
            #https://github.com/3DGenomes/TADbit/blob/ca7cbbb9aa35
    #7591da623bc362cc097fae6de98d/_pytadbit/parsers/hic_bam_parser.py#L278
            strand = 1
            #  Get the data we are interested in
            with open(filepath, 'r') as f:
                header = f.readline()
                for line in f:
                    line = line.split()
                    fragsPos = line[FragsPos].split(';')
                    readId = '-'.join([line[0].split('-')[1], line[0].split('-', 3)[3]])

                    concatemer = []
                    for f in fragsPos:
                        f = f.split(':')
                        chrom = f[0]
                        f2 = f[1].split('-')
                        # in case start is 0
                        start = max(1, int(f2[0]))
                        end = int(f2[1].split('(')[0])
                        seqLen = end - start + 1

                        # save RE fragment
                        concatemer += [[chrCode + str(chrom), str(start),
                               str(strand), str(seqLen), str(start), str(end)]]

                        nConcat += 1

                        if hic_data_style:
                            if removeSingle and len(concatemer) > 1:
                                #  get all posible combinatios between the fragments
                                allComb = list(itertools.combinations(concatemer, 2))
                                nmulti = len(allComb)

                                # Write concatemer interactions
                                for nc, comb in enumerate(allComb, 1):
                                    key = '%s.%s#%s/%s' %(id1, readId,
                                                    nc, nmulti)
                                    toWrite = '%s\t%s\t%s\n' %(key,
                                                            '\t'.join(comb[0]),
                                                            '\t'.join(comb[1]))
                                    # write
                                    with open(outfile, 'a') as f:
                                        f.write(toWrite)
                        else:
                            if removeSingle and len(concatemer) > 1:
                                # add first fragment agains the rest
                                nmulti = len(concatemer) - 1

                                # Write concatemer interactions
                                for nc, comb in enumerate(concatemer[1:], 1):
                                    key = '%s.%s#%s/%s' %(id1, readId,
                                                        nc, nmulti)
                                    toWrite = '%s\t%s\t%s\n' %(key,
                                                            '\t'.join(concatemer[0]),
                                                            '\t'.join(comb))
                                    # write
                                    with open(outfile, 'a') as f:
                                        f.write(toWrite)


            # Store concatemer numbers and concatemers with repeated RE
            nConcatDict[id1]['nConcat'] = nConcat


    return nConcatDict

#  Function to write the TSV files from chiadrop data processed by me
def fromCHIADROPtoTSV_myFormat(files, outPath, all_chromLengths,
                    filterQual=21,
                    hic_data_style=True, chrCode = '',
                    maxFrag=1000000):

    '''
    Function to generate TSV files from the ChIA-drop ones (.multiC)
    :param files: list with paths in which the input ChIA-drop files are located
    :param outPath: String with path in which to store the output tsv files
    :param all_chromLengths: Dictionary with reference genomeas key and
        a text including the lengts of each chromosome like in SAM format
    :param 21 filterQual: Minimum quality accepted
    :param '' chrCode: can be chr for example, in cases in which process data does
        not include it
    :param True hic_data_style: Wether to write the tsv in hic_data style readble
        by TADbit, or just keep minimum posible interactions to be used by our
        pipeline. The first one basically shows all vs all interactions for each
        concatemer, whereas the second one shows the first fragment interaction
        against all the others, saving some space
    :param 1000000 maxFrag: maximum allowed size of fragments to be kept
    '''

    # Set it to de text you want to add before eery chromomse.
    nConcatDict = {}
    #doubleRE = {}
    for filepath in sorted(files):
        cmd = ''
        if filepath.endswith('multiC'):
            id1 = filepath.split('/')[-1][:-7]
            print('## ' + id1 + ' ##')

            if hic_data_style:
                outfile = outPath + '%s_%smaxFr_full.tsv' %(id1, maxFrag)
            else:
                outfile = outPath + '%s_%smaxFr_short.tsv' %(id1, maxFrag)

            #viewPoint = proposedView[id1]

            # get reference genome chromosome coordinates for this sample
            chromLengths = all_chromLengths

            # Create list with concatemers with repeated RE
            #doubleRE[id1] = {}

            # Keep record of number of concatemers check
            nConcat = 0
            nConcatDict[id1] = {}
            # Keep record of number of concatermers with duplicated RE
            #nConcatDup = 0

            ##  First to write, the header of the tsv
            for h in chromLengths.split('\n'):
                if len(h) != 0:
                    cmd += '# CRM %s\n' %'\t'.join(h.split())
            #  Create file and write first lines
            with open(outfile,"w") as f:
                f.write(cmd)

            ##  Then write the data
            # start adding the contacts
            concatemer = []

            # index of positions of interest
            FragsPos = 2
            nFragPos = 1

            # strand position always in forward
            # I believe this should mean forward?
            #https://github.com/3DGenomes/TADbit/blob/ca7cbbb9aa35
    #7591da623bc362cc097fae6de98d/_pytadbit/parsers/hic_bam_parser.py#L278
            strand = 1
            #  Get the data we are interested in
            with open(filepath, 'r') as f:
                #header = f.readline()
                for line in f:
                    line = line.split()
                    nFrag = int(line[nFragPos])
                    if nFrag <= maxFrag:
                        fragsPos = line[FragsPos].split(';')
                        readId = '-'.join(line[0].split('-')[-2:])

                        concatemer = []
                        for f in fragsPos:
                            f = f.split(':')
                            chrom = f[0]
                            f2 = f[1].split('-')
                            # in case start is 0
                            start = max(1, int(f2[0]))
                            end = int(f2[1])
                            seqLen = end - start + 1

                            # save RE fragment
                            concatemer += [[chrCode + str(chrom), str(start),
                                   str(strand), str(seqLen), str(start), str(end)]]

                        nConcat += 1

                        if hic_data_style:
                            if len(concatemer) > 1:
                                #  get all posible combinatios between the fragments
                                allComb = list(itertools.combinations(concatemer, 2))
                                nmulti = len(allComb)

                                # Write concatemer interactions
                                for nc, comb in enumerate(allComb, 1):
                                    key = '%s.%s#%s/%s' %(id1, readId,
                                                    nc, nmulti)
                                    toWrite = '%s\t%s\t%s\n' %(key,
                                                            '\t'.join(comb[0]),
                                                            '\t'.join(comb[1]))
                                    # write
                                    with open(outfile, 'a') as f:
                                        f.write(toWrite)
                        else:
                            if len(concatemer) > 1:
                                # add first fragment agains the rest
                                nmulti = len(concatemer) - 1

                                # Write concatemer interactions
                                for nc, comb in enumerate(concatemer[1:], 1):
                                    key = '%s.%s#%s/%s' %(id1, readId,
                                                        nc, nmulti)
                                    toWrite = '%s\t%s\t%s\n' %(key,
                                                            '\t'.join(concatemer[0]),
                                                            '\t'.join(comb))
                                    # write
                                    with open(outfile, 'a') as f:
                                        f.write(toWrite)


            # Store concatemer numbers and concatemers with repeated RE
            nConcatDict[id1]['nConcat'] = nConcat


    return nConcatDict


# Function to create bedgraphs from the multiC files
def fromCHIADROPtoBedGraph_myFormat(filepath, outPath, all_chromLengths,
                    chrCode = '', maxFrag=1000000, resol=100000, norm=False,
                                   adjustBy=100000000, roundVals=False):

    '''
    Function to generate BedGraph files from the ChIA-drop ones (.multiC)
    :param filepath: path in which the input ChIA-drop file is located
    :param outPath: String with path in which to store the output tsv files
    :param chromList: List with all present chromosomes in the dataset (in
        desired order)
    :param all_chromLengths: dictionary with chromosome names as keys and maximum
        length as integer value
    :param '' chrCode: can be chr for example, in cases in which process data does
        not include it
    :param 1000000 maxFrag: maximum allowed size of fragments to be kept
    :param 100000 resol: resolution at which we will merge reads
    :param False norm: True to normalise by number of total interactions
    :param 100000000 adjustBy: number by which to multiplicate results (if they are too 
        small and you wount to increase the size)
    :param False roundVals: Set to True if you wand to round to 4 decimals. To be
        used only with norm = True 
    '''

    
    def getVal_round(val, normBy):
        return round(val / normBy, 4)
    def getVal_normal(val, normBy):
        return val / normBy
    
    if norm == True and roundVals == True:
        getVal = getVal_round
    else:
        getVal = getVal_normal
        
    # get order of chromosomes
    chromList = sorted(list(all_chromLengths.keys()))
    # Set it to de text you want to add before eery chromomse.
    nConcatDict = {}
    #doubleRE = {}
    allInter = 0
    if filepath.endswith('multiC'):
        id1 = filepath.split('/')[-1][:-7]
        print('## ' + id1 + ' ##')

        outfile = outPath + '%s_%smaxFr.bedgraph' %(id1, maxFrag)
        
        # index of positions of interest
        FragsPos = 2
        nFragPos = 1

        allChromData = {}
        for c in chromList:
            allChromData[c] = defaultdict(int)
        #  Get the data we are interested in binned in specific chunks
        with open(filepath, 'r') as f:
            #header = f.readline()
            for line in f:
                line = line.split()
                nFrag = int(line[nFragPos])
                if nFrag <= maxFrag:
                    fragsPos = line[FragsPos].split(';')

                    for f in fragsPos:
                        f = f.split(':')
                        chrom = f[0]
                        f2 = f[1].split('-')
                        # in case start is 0
                        start = max(1, int(f2[0]))
                        
                        allChromData[chrCode+chrom][start//resol*resol] += 1
                        allInter += 1


        # Now we write the bed file
        if norm == True:
            # we try to multiply positively the result
            normBy = allInter / adjustBy
        else:
            normBy = 1
            
        #print(allInter)
        npeak = 1
        with open(outfile,"w") as fout:
            for chrom in chromList:
                for pos in sorted(allChromData[chrom]):
                    if allChromData[chrom][pos] != 0:
                        # Position is zero based
                        endpos = min(pos + resol, all_chromLengths[chrom])
                        value = getVal(allChromData[chrom][pos], normBy)
                        toWrite = f'{chrom}\t{pos}\t{endpos}\t'
                        toWrite += f'{value}\n'

                        npeak += 1

                        # write
                        fout.write(toWrite)

    return npeak, outfile


## set of functions to write TSV files from our chiadrop data
# Write concatemer interactions
def fullwrite(id1, readId, nmulti, allComb, outfile):
    for nc, comb in enumerate(allComb, 1):
        key = '%s.%s#%s/%s' %(id1, readId, 
                        nc, nmulti)
        toWrite = '%s\t%s\t%s\n' %(key, 
                                '\t'.join(comb[0]), 
                                '\t'.join(comb[1]))
        # write
        with open(outfile, 'a') as f:
            f.write(toWrite) 
        

# Write concatemer interactions
def shortwrite(id1, readId, nmulti, concatemer, outfile):
    for nc, comb in enumerate(concatemer[1:], 1):
        key = '%s.%s#%s/%s' %(id1, readId,
                            nc, nmulti)
        toWrite = '%s\t%s\t%s\n' %(key,
                                '\t'.join(concatemer[0]),
                                '\t'.join(comb))
        # write
        with open(outfile, 'a') as f:
            f.write(toWrite)
            
def writeGEM(hic_data_style, concatemer, id1, prevGEMid,
                            outfile):
    ## Write
    if hic_data_style:
        #  get all posible combinatios between the fragments
        allComb = list(itertools.combinations(concatemer, 2))
        nmulti = len(allComb)

        # Write concatemer interactions
        fullwrite(id1, prevGEMid, nmulti, allComb, outfile)

    else:
        # add first fragment agains the rest
        nmulti = len(concatemer) - 1

        # Write concatemer interactions
        shortwrite(id1, prevGEMid, nmulti, concatemer, outfile)

#  Function to write the TSV files from chiadrop files
def fromCHIADROPtoTSV(outPath, all_chromLengths,
                    filterQual=21, 
                    hic_data_style=True, chrCode = ''):

    '''
    Function to generate TSV files from the HDF5 ones
    param outPath: String with path in which to store the output tsv files
    param all_chromLengths: Dictionary with reference genomeas key and
        a text including the lengts of each chromosome like in SAM format
    param 21 filterQual: Minimum quality accepted
    :param '' chrCode: can be chr for example, in cases in which process data does 
        not include it
    param True hic_data_style: Wether to write the tsv in hic_data style readble
        by TADbit, or just keep minimum posible interactions to be used by our 
        pipeline. The first one basically shows all vs all interactions for each
        concatemer, whereas the second one shows the first fragment interaction
        against all the others, saving some space
    '''
    
    files = os.listdir(outPath)
    # Set it to de text you want to add before eery chromomse.
    nConcatDict = {}
    #doubleRE = {}
    for filepath in sorted(files):
        cmd = ''
        if filepath.endswith('region'):
            id1 = '_'.join(filepath.split('_',2)[0:2])
            print('## ' + id1 + ' ##')
            if not os.path.exists(outPath + 'tsv/'):
                os.makedirs(outPath + 'tsv/')

            if hic_data_style:
                outfile = outPath + 'tsv/%s_full.tsv' %(id1)
            else:
                outfile = outPath + 'tsv/%s_short.tsv' %(id1)

            #viewPoint = proposedView[id1]

            # get reference genome chromosome coordinates for this sample
            chromLengths = all_chromLengths

            # Create list with concatemers with repeated RE
            #doubleRE[id1] = {}

            # Keep record of number of concatemers check
            nConcat = 0
            nConcatDict[id1] = {}
            # Keep record of number of concatermers with duplicated RE
            #nConcatDup = 0

            ##  First to write, the header of the tsv
            for h in chromLengths.split('\n'):
                if len(h) != 0:
                    cmd += '# CRM %s\n' %'\t'.join(h.split())
            #  Create file and write first lines
            with open(outfile,"w") as f:
                f.write(cmd)

            ##  Then write the data
            # start adding the contacts
            concatemer = []

            # index of positions of interest
            nFragPos = 3
            chromPos = 0
            startPos = 1
            endPos = 2
            GEMidPos = 4


            # get first GEM id
            with open(outPath + filepath, 'r') as f:
                first_line = f.readline()
                prevGEMid = first_line.split()[GEMidPos]
                nFrag = int(first_line.split()[nFragPos])

            # strand position always in forward
            # I believe this should mean forward?
            #https://github.com/3DGenomes/TADbit/blob/ca7cbbb9aa35
    #7591da623bc362cc097fae6de98d/_pytadbit/parsers/hic_bam_parser.py#L278
            strand = 1
            #  Get the data we are interested in
            with open(outPath + filepath, 'r') as f:
                #header = f.readline()
                for line in f:
                    line = line.split()

                    GEMid = line[GEMidPos]

                    # if we are in a new concatemer
                    if prevGEMid != GEMid:
                        writeGEM(hic_data_style, concatemer, id1, prevGEMid,
                                outfile)

                        ## restart parameters 
                        # restart list with GEM fragments coordiantes
                        concatemer = []
                        # reset previous id to actual one
                        prevGEMid = GEMid
                        nConcat += 1

                    # get coordinates and info
                    chrom = line[chromPos]
                    start = int(line[startPos])
                    end = int(line[endPos])
                    nFrag = int(line[nFragPos])
                    seqLen = end - start + 1

                    # save GEM
                    concatemer += [[chrCode + str(chrom), str(start),
                           str(strand), str(seqLen), str(start), str(end)]]

                # write last GEM
                writeGEM(hic_data_style, concatemer, id1, prevGEMid,
                            outfile)

                ## restart parameters 
                # restart list with GEM fragments coordiantes
                concatemer = []
                # reset previous id to actual one
                prevGEMid = GEMid
                nConcat += 1

            # Store concatemer numbers and concatemers with repeated RE
            nConcatDict[id1]['nConcat'] = nConcat


    return nConcatDict


######  FUNCTIONS FOR THE NORMALISATON  ######

# Obtain multiContacts from file just in a pairwise manner
def goThroughReads_yesDiag(section_pos, line, interPerBin,
                                resol, nRead=0):
    '''
    Function to obtain interaction frecuencies per bin from 
        normal TSV with no multiContact data spected
        
    :param section_pos: Dictionary with the chromosomes
        reference name as keys and a list or tuple of
        the range of bins in which they lie. Ej.:
        {'chr1':(0, 5000)}. Is 0 index, and last bin 
        from range, is not included inside, so next
        chromosome could be {'chr2':(5000, 10000)}
    :param line: list or tuple with:
        [chr1, startPos1, chr2, startPos2]
    :param interPerBin: defaultdict(int) with the number
        of times each bin interacts
    :param resol: Resolution at wich we are going to be 
        normalising our data
    :param 0 nRead: Integer indicating number of reads 
        counted
    
    '''
   
    # store each apparition of a fragment in a concatemer
    #we store the bin of the mapping start position
    fragment1 = (int(line[1]) // resol) + section_pos[line[0]][0]
    interPerBin[fragment1] += 1
    
    fragment2 = (int(line[3]) // resol) + section_pos[line[2]][0]
    interPerBin[fragment2] += 1
    
    # New concatemer seen
    nRead += 1
    
    return interPerBin, nRead


# Obtain multiContacts from file just in a pairwise manner. Remove diagonal
def goThroughReads_noDiag(section_pos, line, interPerBin,
                                resol, nRead=0):
    '''
    Function to obtain interaction frecuencies per bin from 
        normal TSV with no multiContact data spected
        
    :param section_pos: Dictionary with the chromosomes
        reference name as keys and a list or tuple of
        the range of bins in which they lie. Ej.:
        {'chr1':(0, 5000)}. Is 0 index, and last bin 
        from range, is not included inside, so next
        chromosome could be {'chr2':(5000, 10000)}
    :param line: list or tuple with:
        [chr1, startPos1, chr2, startPos2]
    :param interPerBin: defaultdict(int) with the number
        of times each bin interacts
    :param resol: Resolution at wich we are going to be 
        normalising our data
    :param 0 nRead: Integer indicating number of reads 
        counted
    
    '''
   
    # store each apparition of a fragment in a concatemer
    # we store the bin of the mapping start position
    fragment1 = (int(line[1]) // resol) + section_pos[line[0]][0]
    fragment2 = (int(line[3]) // resol) + section_pos[line[2]][0]
    
    if fragment1 != fragment2:
        interPerBin[fragment1] += 1
        interPerBin[fragment2] += 1
        # New concatemer seen
        nRead += 1
    
    return interPerBin, nRead

# function to get the bins list per chromosome
def getSectionPos(infile, resol):
    bamfile = AlignmentFile(infile, 'rb')
    bam_refs = bamfile.references
    bam_lengths = bamfile.lengths

    sections = OrderedDict(list(zip(bam_refs,
                               [x for x in bam_lengths])))
    total = 0
    section_pos = OrderedDict()
    for crm in sections:
        section_pos[crm] = (total, total + (sections[crm] // resol + 1))
        total += (sections[crm] // resol + 1)
        
    bamfile.close()
    return section_pos


# Open tsv and obtain contact frecuencies in a pairwise manner
def getInteractionsPerBin(infile, resol, locusCh=False,
                         regRange = False, returnNread = False,
                         filter_exclude=(1, 2, 3, 4, 6, 7, 8, 9, 10),
                         diagonal=True):
    '''
    Function to get the number of concatemers were a bin of interest is 
        appearing (is bin based, so all fragments which start inside of
        a bin margin will be joined)
        
    :param infile: Path to the input file. If TADbit style TSV, be sure 
        it ends with .tsv. If usual BAM file, be sure it is sorted, the
        index is located in the same folder, and it ends with .bam
    :param resol: Resolution at wich we are going to be normalising our
        data
    :param False locusCh: Set to True if you want to return just data to
        normalise the bin between the smallest and biggest binned coordinate
        in regRange. Even if True the whole bam file must be cheked, so
        wont safe any time.
    :param False regRange: list or tuple with the first and last binned
        coordinates we wont to retrieve.
    :param False returnNread: wether you want or not nRead to be
        returned. It returns the number of reads check.
    :param (1, 2, 3, 4, 6, 7, 8, 9, 10) filter exclude: filters to define the
        set of valid pair of reads. Just valid for TADbit style BAMfiles. If
        working with already filtered non TADbit BAM set filter_exclude = ()
    :param True diagonal: True if you want to count diagonal interactions
    '''
    
    # change filter exclude to binary
    if not isinstance(filter_exclude, int):
        filter_exclude = filters_to_bin(filter_exclude)

    if diagonal:
        goThroughReads = goThroughReads_yesDiag
    else:
        goThroughReads = goThroughReads_noDiag

    interPerBin = defaultdict(int)
    
    nRead = 0
    
    if infile.endswith('.bam'):
        # get section positions 
        section_pos = getSectionPos(infile, resol)

        # check if this BAM file is not TADbit style
        bamfile = AlignmentFile(infile, 'rb')
        if 'Hicup' in bamfile.text:
            print('It seems this BAM file was produced outside TADbit, make \
sure it has already been filtered')
            if filter_exclude != ():
                print('Consider changing filter_exclude so its value is () \
or you might get no reads')
                
        bamfile.close()

        bamfile = AlignmentFile(infile, 'rb')
        for r in bamfile.fetch():
            # Check if it is among positions to be filtered
            # BEWARE that it follows TADbit format
            if r.flag & filter_exclude:
                continue
            crm1 = r.reference_name
            pos1 = r.reference_start + 1
            crm2 = r.next_reference_name
            pos2 = r.next_reference_start + 1
            
            line = [crm1, pos1, crm2, pos2]
            interPerBin, nRead = goThroughReads(section_pos, line, interPerBin,
                                    resol, nRead=nRead)
        bamfile.close()
    

    elif infile.endswith('.tsv'):
        with open(infile, 'r') as f:
            # get chromosome lengths
            sections = OrderedDict()
            while True:
                line = f.readline()
                if line.startswith('#'):
                    line = line.split()
                    if line[1] == 'CRM':
                        sections[line[2]] = int(line[3])
                        
                elif line.startswith('@'):
                    pass
                else:
                    break
                    
            # create binned positioning for chromosomes
            total = 0
            section_pos = OrderedDict()
            for crm in sections:
                section_pos[crm] = (total, total + (sections[crm] // resol + 1))
                total += (sections[crm] // resol + 1)

            # iterate over reads
            line = line.split()
            line = [line[1], line[2], line[7], line[8]]

            # Run current line (first one)
            interPerBin, nRead = goThroughReads(section_pos, line, interPerBin,
                                    resol, nRead=nRead)

            # Go for next lines
            for line in f:
                line = line.split()
                line = [line[1], line[2], line[7], line[8]]
                
                interPerBin, nRead = goThroughReads(section_pos, line, interPerBin,
                                    resol, nRead=nRead)

    
     # Get genomic coordinates for region of interest
    if locusCh == True:
        regionStart = min(regRange)
        regionEnd = max(regRange) 

        ## modify if we want data from all genome
        # Get all the fragments that start inside this coordinates
        keys = [k for k in interPerBin.keys() if (regionStart <= 
                                                            k <= 
                                                            regionEnd)]
        regInterPerBin = defaultdict(int)
        for k in keys:
            regInterPerBin[k] += interPerBin[k]
            
    # Or not
    else:
        regInterPerBin = interPerBin
            

    if returnNread == False:
        return regInterPerBin
    else:
        return regInterPerBin, nRead


# Normalise by frecuencies given the presence of each interacting fragment                                           
def frecuenciesNorm(hic_data, resol, regRange, concatemersBin, multResult=100, 
                    keep=False, mininter=0, positionAdjust=0):
    
    '''
    param 0 positionAdjust: In case the positions from concatemersBin are taking 
        into account bining from the whole genome, but we just load in hic_data
        one chromosome. Here concatemersBin will be substracted added to the 
        positions in regRange in order to compensate this
    
    '''
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
                    #if concatemersBin[bin1] == 0:
                    #    if concatemersBin[bin2] == 0:
                    #        divider = 1
                    #    else:
                    #        divider = concatemersBin[bin2]
                    #elif concatemersBin[bin2] == 0:
                    #    divider = concatemersBin[bin1]
                    #else:
                    divider = (concatemersBin[bin1 + positionAdjust] + 
                               concatemersBin[bin2 + positionAdjust] - 
                               hic_data[bin1, bin2])

                    #divider = float(concatemersBin[bin1] + concatemersBin[bin2])

                    #if divider == 0:
                    #    divider = 1
                    # if divider is 0 it means concatemersBin was taken with another index
                    #ie just checking a chromosome, or whole genome and here just 
                    #normalising a file were we loaded one chromosome
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
                #if concatemersBin[ke[0]] == 0:
                #    if concatemersBin[ke[1]] == 0:
                #        divider = 1
                #    else:
                #        divider = concatemersBin[ke[1]]
                #elif concatemersBin[ke[1]] == 0:
                #    divider = concatemersBin[ke[0]]
                #else:
                divider = (concatemersBin[ke[0] + positionAdjust] + 
                               concatemersBin[ke[1] + positionAdjust] - 
                               hic_data[ke[0], ke[1]])

                #divider = float(concatemersBin[bin1] + concatemersBin[bin2])

                #if divider == 0:
                #    divider = 1
                # if both are zero 
                norm_data[ke[0], ke[1]] = hic_data[ke[0], ke[1]] / float(divider)
                norm_data[ke[1], ke[0]] = norm_data[ke[0], ke[1]]
        
    return norm_data