import argparse
from taddyn.modelling.impoptimizer  import IMPoptimizer
from taddyn import Chromosome, Experiment
import numpy as np
#import seaborn as sns
import sys
import os, errno
import shutil
import random

def PCHiC_filtering(exp, index=0):
    zeros = {}
    for i in range(exp.size):
        lineAdd = 0
        for j in xrange(exp.size):
            lineAdd += exp.norm[index]['matrix'].get(i * exp.size + j, 0)
        if lineAdd == 0:
            zeros[i] = 0
    exp._zeros = [zeros]

def cToDot(text):
    newText = ''
    for t in text:
        if t == 'c':
            newText += '.'
        elif t == 'n':
            newText += '-'
        else:
            newText += t
    return newText


parser = argparse.ArgumentParser(description='')
parser.add_argument('-l','--lowfreq', help='lowfreq',required=True)
parser.add_argument('-m','--maxdist',help='maxdist_range', required=True)
parser.add_argument('-d','--dcutoff_range',help='dcutoff_range', required=True)
parser.add_argument('-u','--upperfreq',help='upperfreq', required=True)
parser.add_argument('-p','--pathtomtrx',help='path_to_matrix', required=True)
parser.add_argument('-t','--jobtime',help='jobtime_HH:MM:SS', required=True)
parser.add_argument('-nm','--nmodels',help='nmodels', required=True)
parser.add_argument('-tp','--temp_path',help='path_to_tmp_files', required=False)

args = parser.parse_args()
lowfreq= float(args.lowfreq)
uperfreq= float(args.upperfreq)
try:
    maxdist= int(args.maxdist)
except:
    maxdist= args.maxdist
    print 'We dont consider floating points in maxdist: %s' %(maxdist)
    maxdist= int(float(maxdist))

dcutoff_range= [float(s) for s in cToDot(args.dcutoff_range).split('_')]
matPath=args.pathtomtrx
jobTime=args.jobtime
nmodels = int(args.nmodels)
tempOut=args.temp_path
#matPath='.../MatrixFreqNorm_WPL-KOD_reg4_chr8-120780000-122030000_10000kb'

flag = matPath.split('/')[-2]
path='/'.join(matPath.split('/')[:-1]) + '/'
res = int(matPath.split('_')[-1][:-2])

# select dcutofs
dcutoff_range=np.arange(dcutoff_range[0],dcutoff_range[1],dcutoff_range[2]) #Cutoff
dcutoff_range2 = []
for dc in dcutoff_range:
    if maxdist - dc >= -50:
        dcutoff_range2.append(float(dc))

#Rao_HIC/load8/Matrix_28331000_28568000
test_chr = Chromosome(name='Test',centromere_search=False,
                      species='Homo sapiens', assembly='37')#, max_tad_size=260000)
test_chr.add_experiment('test',exp_type='Hi-C', identifier='GM128Rao', resolution=res,
                        norm_data=matPath)


exp = test_chr.experiments[0]
print exp
print exp.size


## Filtering
# If HiC data
if 'OneD' in matPath:
    exp.filter_columns(silent=False,draw_hist=False)
# If pcHiC or virtual pcHiC
else:
    PCHiC_filtering(exp)


##############################
keep_restart_out_dir = path + 'lammpsSteps/' 
if not os.path.exists(keep_restart_out_dir):
    os.makedirs(keep_restart_out_dir)
jobName = 'LF%sUF%sMdis%s_%sbp'%(str(lowfreq),str(uperfreq),str(maxdist), str(res))
keep_restart_out_dir = path + 'lammpsSteps/jobArray_%s/' %jobName
if not os.path.exists(keep_restart_out_dir):
        os.makedirs(keep_restart_out_dir)
restart_path = keep_restart_out_dir

#keep_restart_out_dir = None
#restart_path = False

# define initial seed in order it gets totally different or same models
#initial_seed = random.choice(range(0, 100000000, nmodels))
# cant use random because then the folder name to look for the individual
#models changes lammps_n
initial_seed = 0

dcut_text = '-'.join(str(d) for d in dcutoff_range)
optimizer = IMPoptimizer(exp, start=1, end=exp.size, n_models=nmodels, n_keep=nmodels,  tool='lammps', tmp_folder= tempOut)
optimizer.run_grid_search(n_cpus=min(nmodels, 8), lowfreq_range=[float(lowfreq)],
                          upfreq_range=[float(uperfreq)],
                          maxdist_range=[float(maxdist)],
                          dcutoff_range=dcutoff_range2,
                          scale_range=[0.01][:], verbose=True,
                          timeout_job=jobTime,
			  #savedata=path+'opt_LF%sUF%sMdis%s_%sbp.modelsPickle'%(str(lowfreq),str(uperfreq),str(maxdist), str(res)),
			  cleanup=True, hide_log=True,
              initial_seed=initial_seed,
              initial_conformation='random',
              # Lines to make timePoints and load them if they exist
              keep_restart_out_dir=keep_restart_out_dir,
              restart_path=restart_path,
              store_n_steps=10)

outfile=path+'opt_LF%sUF%sC%sMdis%s_%sbp.txt'%(str(lowfreq),str(uperfreq),dcut_text,str(maxdist), str(res))
if keep_restart_out_dir != None:
    # Now we need to check if all models finished before storing any output
    modelFiles = os.listdir(keep_restart_out_dir)
    nmodelFound = 0
    for mo in modelFiles:
        k = mo.split('_')[1]
        if os.path.exists('%s/%s/finishedModel_%s.pickle' %(keep_restart_out_dir, mo, k)):
            nmodelFound += 1
    if nmodelFound == nmodels:
        print 'All models finished correctly'
        optimizer.write_result(outfile)
        # remove steps to save disk quota
        shutil.rmtree(keep_restart_out_dir) 
    else:
        print 'Some models didnt finished'
else:
    print 'Beware that some models might have failed, check WARNING prints'
    optimizer.write_result(outfile)

