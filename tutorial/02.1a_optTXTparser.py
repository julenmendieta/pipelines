import sys
import matplotlib
#matplotlib.use('Agg')
import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse
from mpl_toolkits.mplot3d import Axes3D

parser = argparse.ArgumentParser(description='')
parser.add_argument('-p','--pathOut',help='path_to_output', required=True)
parser.add_argument('-j','--joinWithOld',help='join_with_old_results', required=True)

args = parser.parse_args()
path=args.pathOut
join=args.joinWithOld

print '#' * 30
print path

# get data from previous files
oldFiles = {}
for file1 in os.listdir(path):
    if file1.startswith('Val') and file1.endswith('.txt'):
        dcut = 'C' + str(float(file1.split('C')[-1][:-4]))
        oldFiles[dcut] = []
        with open(path + file1, 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    oldFiles[dcut].append(line)


#path = './'
header=''
val_dict={}
for file in os.listdir(path):
    if file.startswith('opt') and file.endswith('.txt'):
        fclean=file[4:-4].replace('LF',':LF').replace('UF',':UF').replace('C',':C').replace('Mdis',':Mdis')
        #opt_LF-0.9UF0.1C2Mdis1400.txt
        c=fclean.split(':')[3]#cutoff
        firstC = 'C' + str(float(c.split('-')[0].split('C')[1]))
        allC = [firstC] + ['C' + str(float(cc)) for cc in c.split('-')[1:]]
        for c1 in allC:
            val_dict[c1]=[]

for file in os.listdir(path):
    if file.startswith('opt') and file.endswith('.txt'):
            fop=open(path + file,'r')
            fclean=file[4:-4].replace('LF',':LF').replace('UF',':UF').replace('C',':C').replace('Mdis',':Mdis')
            for l in fop.readlines():
                if l.startswith('##'):
                    header=l
                else:
                    if l.startswith('#'):
                        header += l
                    else:
                        # remove files with zero values
                        corr = l.split('\t')[-1]
                        if corr != corr or float(corr) == 0:
                            print 'BAD: %s' %file
                        else:
                            val_dict['C' + str(float(l.split()[5]))].append(l)

#print header

# move to oldsFile or join it 
if join == 'yes':
    for dcut in oldFiles.keys():
        # if we had check again this distance cutoff
        if dcut in val_dict:
            for line1 in oldFiles[dcut]:
                present = False
                for line2 in val_dict[dcut]:
                    if line1.split()[:-1] == line2.split()[:-1]:
                        present = True
                # if we dont have this region, we add it
                if present == False:
                    val_dict[dcut].append(line1)
        # if we dont have any data of this dcut, we store it full
        else:
            val_dict[dcut] = oldFiles[dcut]
else:
    # store it separated otherwise
    for v in oldFiles.keys():
            fileout=open("%sValOptimisation_old_%s.txt"%(path, v),'w')
            fileout.write(header)
            print v,len(oldFiles[v])
            for ll in oldFiles[v]:
                fileout.write(ll)
            fileout.close()

# store new data
for v in val_dict.keys():
    fileout=open("%sValOptimisation%s.txt"%(path, v),'w')
    fileout.write(header)
    print v,len(val_dict[v])
    for ll in val_dict[v]:
        fileout.write(ll)
    fileout.close()

