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

args = parser.parse_args()
path=args.pathOut

#path = './'
header=''
val_dict={}
for file in os.listdir(path):
	if file.startswith('opt') and file.endswith('.txt'):
			fclean=file[4:-4].replace('LF',':LF').replace('UF',':UF').replace('C',':C').replace('Mdis',':Mdis')
			#opt_LF-0.9UF0.1C2Mdis1400.txt
			c=fclean.split(':')[3]#cutoff
			val_dict[c]=[]

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
						val_dict[fclean.split(':')[3]].append(l)

print header


for v in val_dict.keys():
	fileout=open("%sValOptimisation%s.txt"%(path, v),'w')
	fileout.write(header)
	print v,len(val_dict[v])
	for ll in val_dict[v]:
		fileout.write(ll)
	fileout.close()

