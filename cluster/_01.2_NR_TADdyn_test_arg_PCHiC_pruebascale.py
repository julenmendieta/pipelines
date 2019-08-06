import argparse
from pytadbit.modelling.impoptimizer  import IMPoptimizer
from pytadbit import Chromosome, Experiment
import numpy as np
#import seaborn as sns
import sys
sys.path.append('/home/devel/jmendietaesteban/PCHIC/codigo')
sys.path.append('/home/devel/jmendietaesteban/programas/PhD/3C')
import handy3C
import PCHiC
import matplotlib.pyplot as plt
import os, errno
import shutil
from pebble import ProcessPool

parser = argparse.ArgumentParser(description='')
parser.add_argument('-l','--lowfreq', help='lowfreq',required=True)
parser.add_argument('-m','--maxdist',help='maxdist_range', required=True)
parser.add_argument('-d','--dcutoff_range',help='dcutoff_range', required=True)
parser.add_argument('-u','--upperfreq',help='upperfreq', required=True)
parser.add_argument('-p','--pathtomtrx',help='path_to_matrix', required=True)
parser.add_argument('-t','--jobtime',help='jobtime_HH:MM:SS', required=True)
parser.add_argument('-nm','--nmodels',help='nmodels', required=True)

args = parser.parse_args()
lowfreq= args.lowfreq
uperfreq= args.upperfreq
maxdist= args.maxdist
dcutoff_range= args.dcutoff_range
matPath=args.pathtomtrx
jobTime=args.jobtime
nmodels = int(args.nmodels)
#matPath='/scratch/devel/jmendieta/deLaat/modelling/reg4_WPL-KOD/MatrixFreqNorm_WPL-KOD_reg4_chr8-120780000-122030000_10000kb'

flag = matPath.split('/')[-2]
tempOut = "/scratch_tmp/%s/temp" %(flag)
path='/'.join(matPath.split('/')[:-1]) + '/'
res = int(matPath.split('_')[-1][:-2])


#Rao_HIC/load8/Matrix_28331000_28568000
test_chr = Chromosome(name='Test',centromere_search=False,
                      species='Homo sapiens', assembly='37')#, max_tad_size=260000)
# If interaction index and not matrix
if matPath.split('/')[-1].startswith('interaction'):
	regionStart = int(matPath.split('/')[-1].split('_')[3].split('-')[1])
	regionEnd = int(matPath.split('/')[-1].split('_')[3].split('-')[2])
	matPath = handy3C.tableToMatrix(matPath, regionStart, regionEnd, resol=res, dirout="")
	test_chr.add_experiment('test',exp_type='Hi-C', identifier='GM128Rao', resolution=res,
                        norm_data=[matPath])
else:
	test_chr.add_experiment('test',exp_type='Hi-C', identifier='GM128Rao', resolution=res,
        	                norm_data=matPath)


print '####'


exp = test_chr.experiments[0]
print exp
print exp.size




#print test_chr._find_centromere(exp)

#exp.view(normalized=False)#,where = 'down')

PCHiC.PCHiC_filtering(exp)
#exp.filter_columns(silent=False,draw_hist=False, perc_zero=100)
#exp.normalize_hic(silent=False, factor=None)

#fig=plt.figures()
#sns.distplot(z_score)
#plt.savefig(flag+'_Zscore.pdf')

# I want the modelling to stop in the time limit i said
time1 = [int(t) for t in jobTime.split(':')]
time1 = (time1[0] * 3600) + (time1[1] * 60) + time1[2]

optimizer = IMPoptimizer(exp, start=1, end=exp.size, n_models=nmodels, n_keep=nmodels,  tool='lammps', tmp_folder= tempOut)
optimizer.run_grid_search(n_cpus=min(nmodels, 8), lowfreq_range=[float(lowfreq)],
                          upfreq_range=[float(uperfreq)],
                          maxdist_range=[float(maxdist)],
                          dcutoff_range=[float(dcutoff_range)],
                          scale_range=[0.01][:], verbose=3,
                          timeout_job=time1,
			  savedata=path+'opt_LF%sUF%sC%sMdis%s_%sbp.models'%(str(lowfreq),str(uperfreq),str(dcutoff_range),str(maxdist), str(res)),
			  kfactor=0.1, cleanup=True)

outfile=path+'opt_LF%sUF%sC%sMdis%s_%sbp.txt'%(str(lowfreq),str(uperfreq),str(dcutoff_range),str(maxdist), str(res))
optimizer.write_result(outfile)

