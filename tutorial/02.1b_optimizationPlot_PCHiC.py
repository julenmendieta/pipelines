from taddyn.modelling.impoptimizer  import IMPoptimizer
from taddyn import Chromosome
import os
import argparse
import parser

def PCHiC_filtering(exp, index=0):
    zeros = {}
    for i in range(exp.size):
        lineAdd = 0
        for j in xrange(exp.size):
            lineAdd += exp.norm[index]['matrix'].get(i * exp.size + j, 0)
        if lineAdd == 0:
            zeros[i] = 0
    exp._zeros = [zeros]

parser = argparse.ArgumentParser(description='')
parser.add_argument('-p','--pathtomtrx',help='path_to_matrix', required=True)
parser.add_argument('-nm','--nmodels',help='nmodels', required=False)
parser.add_argument('-nkeep','--nkeep',help='keep_n_models', required=False)

args = parser.parse_args()
matrixPath=args.pathtomtrx
nmodels = args.nmodels
nkeep = args.nkeep


#matrixPath = "/scratch/devel/jmendieta/Ferrer/modelling/PIdyn/reg32Norm/MatrixNormMin_43535000_44430000_compPI"
#outpath = "/scratch/devel/jmendieta/Ferrer/modelling/PIdyn/reg32Norm/scale01/" # can be ""

chrin = ''
outpath='/'.join(matrixPath.split('/')[:-1]) + '/'
res = int(matrixPath.split('_')[-1][:-2])
if nmodels == None:
    nmodels = 100
if nkeep == None:
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
test_chr.add_experiment('test',exp_type='Hi-C', resolution=res,
                         norm_data=matrixPath)


exp = test_chr.experiments[0]

# If HiC data
if 'OneD' in matrixPath:
        exp.filter_columns(silent=False,draw_hist=False)
# If pcHiC or virtual pcHiC
else:
        PCHiC_filtering(exp)


for opt in optFpath:
	print opt
	optim=IMPoptimizer(exp,start=1, end=exp.size, close_bins=1, n_models=nmodels, n_keep=nkeep)
	optim.load_from_file(opt)
	optim.plot_2d(savefig=opt[:-3]+'pdf',show_best=10)#[0.2,0.8])#skip={'maxdist':2000}
	print optim.get_best_parameters_dict()

