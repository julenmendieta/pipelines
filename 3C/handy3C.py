import pandas as pd


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
