import argparse
#from pytadbit.modelling.impoptimizer import IMPoptimizer
from pytadbit.modelling.impoptimizer  import IMPoptimizer
from pytadbit import Chromosome
import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys
# CRG
sys.path.append('/users/mmarti/jmendieta/codigo')
sys.path.append('/users/mmarti/jmendieta/programas/PhD/3C')
# CNAG
sys.path.append('/home/devel/jmendietaesteban/PCHIC/codigo')
sys.path.append('/home/devel/jmendietaesteban/programas/PhD/3C')
import handy3C
import PCHiC
import glob
import os
import parser

parser = argparse.ArgumentParser(description='')
parser.add_argument('-p','--pathtomtrx',help='path_to_matrix', required=True)

args = parser.parse_args()
matrixPath=args.pathtomtrx

#matrixPath = "/scratch/devel/jmendieta/Ferrer/modelling/PIdyn/reg32Norm/MatrixNormMin_43535000_44430000_compPI"
#outpath = "/scratch/devel/jmendieta/Ferrer/modelling/PIdyn/reg32Norm/scale01/" # can be ""

chrin = ''
outpath='/'.join(matrixPath.split('/')[:-1]) + '/'
res = int(matrixPath.split('_')[-1][:-2])
nmodels = 100
nkeep = 100
#res = 5000


optFpath = []
optFiles = os.listdir(outpath)
for opt in optFiles:
	if opt.startswith('ValOptimisation') and not opt.endswith('pdf'):
		optFpath.append(outpath + opt)


test_chr = Chromosome(name='Test%s'%chrin,centromere_search=False,
                      species='Homo sapiens', assembly='na')#, max_tad_size=260000)
# If interaction index and not matrix
if matrixPath.split('/')[-1].startswith('interaction'):
        regionStart = int(matrixPath.split('/')[-1].split('_')[3].split('-')[1])
        regionEnd = int(matrixPath.split('/')[-1].split('_')[3].split('-')[2])
        matrixPath = handy3C.tableToMatrix(matrixPath, regionStart, regionEnd, resol=res, dirout="")
        test_chr.add_experiment('test',exp_type='Hi-C', resolution=res,
                        norm_data=[matrixPath])
else:
        test_chr.add_experiment('test',exp_type='Hi-C', resolution=res,
                                norm_data=matrixPath)


exp = test_chr.experiments[0]

# If HiC data
if 'OneD' in matrixPath:
        exp.filter_columns(silent=False,draw_hist=False)
# If pcHiC or virtual pcHiC
else:
        PCHiC.PCHiC_filtering(exp)


for opt in optFpath:
	print opt
	optim=IMPoptimizer(exp,start=1, end=exp.size, close_bins=1, n_models=nmodels, n_keep=nkeep)
	optim.load_from_file(opt)
	optim.plot_2d(savefig=opt[:-3]+'pdf',show_best=10)#[0.2,0.8])#skip={'maxdist':2000}
	print optim.get_best_parameters_dict()

	# Ya que estamos mostramos los modelos
	# dic = optim.get_best_parameters_dict()
	# file1 = 'opt_LF' + str(dic['lowfreq']) + 'UF' + str(dic['upfreq']) + 'C2Mdis' + str(int(dic['maxdist'])) + '.txt'
	# #print file1
	# if dic['scale'] == 0.01:
	# 	files2 = glob.glob('*Mdis100.txt')
	# elif dic['scale'] == 0.005:
	# 	files2 = glob.glob('*Mdis50.txt')
	# best = ""
	# bestCor = 0
	# for fi in files2:
	# 	with open(fi, 'r') as f:
	# 		for i, line in enumerate(f):
	# 			# In second line take header
	# 			if i == 1:
	# 				line = line.split('\t')
	# 				posLF = line.index('low_freq')
	# 				posUF = line.index('up_freq')
	# 			# Compare values
	# 			elif i == 2:
	# 				line = line.split('\t')
	# 				if line[5] > bestCor and (line[posLF] <= line[posUF]):
	# 					bestCor = line[5][:-1]
	# 					best = fi
	# print best

