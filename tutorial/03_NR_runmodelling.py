from taddyn import Chromosome
import sys
import argparse
import numpy as np
#import seaborn as sns
import os, errno
import shutil
import random
from taddyn.utils.modelAnalysis     import save_models


#import argparse
# Define function to create lammps output folder
def mkdir(dnam):
    try:
        os.mkdir(dnam)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(dnam):
            pass
        else:
            raise

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
parser.add_argument('-l','--lowfreq', help='lowfreq',required=True)
parser.add_argument('-m','--maxdist',help='maxdist_range', required=True)
parser.add_argument('-d','--dcutoff_range',help='dcutoff_range', required=True)
parser.add_argument('-u','--upperfreq',help='upperfreq', required=True)
parser.add_argument('-lf','--lammpsfolder',help='folder_temp_lammps', required=True)
parser.add_argument('-p','--pathtomtrx',help='path_to_matrix', required=True)
parser.add_argument('-t','--jobtime',help='jobtime_HH:MM:SS', required=True)
parser.add_argument('-tp','--temp_path',help='path_to_tmp_files', required=False)
parser.add_argument('-nm','--nmodels',help='model_to_build', required=False)
parser.add_argument('-ncpu','--ncpu',help='number_cpus', required=False)
parser.add_argument('-po','--pathOut',help='output_models_path', required=False)



args = parser.parse_args()
low=float(args.lowfreq)
up=float(args.upperfreq)
try:
    maxd= int(args.maxdist)
except:
    maxd= args.maxdist
    print 'We dont consider floating points in maxdist: %s' %(maxd)
    maxd= int(float(maxd))
c=float(args.dcutoff_range)
lammpsOut=args.lammpsfolder
matPath=args.pathtomtrx
jobTime=args.jobtime
tempOut=args.temp_path
nmodels=args.nmodels
n_cpu=args.ncpu
pathOut=args.pathOut

res = int(matPath.split('_')[-1][:-2])
lampsFlag = lammpsOut.split('_')[-1]
# need to mantain this for the clustering scripts
try:
    a = int(lampsFlag)
except:
    lampsFlag = 1

chrin='na'
path='/'.join(matPath.split('/')[:-1]) + '/'
res = int(matPath.split('_')[-1][:-2])
cell = matPath.split('/')[-1].split('_')[1]
region = matPath.split('/')[-1].split('_')[2]
flag = '%s_%s' %(cell, region)

flag_name='%s_C%sL%sU%sM%sRes%s' %(flag, c, low, up, maxd, res)

# if we provide an alternative output dir we change it 
if tempOut != None:
    lammpsOut = tempOut

if n_cpu == None:
    n_cpu = 1
else:
    n_cpu = int(n_cpu)
# check if we asked for more than one model
if nmodels == None:
    nmodels = 1
else:
    nmodels = int(nmodels)
if pathOut == None:
    pathOut= path + 'finalModel/'

# Move to lammpsOut directory so its output gets stored there
#mkdir(lammpsOut)
#os.chdir(lammpsOut)


#matPath='/scratch/devel/jmendieta/Ferrer/modelling/PIdyn/reg32Norm/MatrixNormMin_43535000_44430000_compPI'

# Create outdir if does not exist
if not os.path.exists(pathOut):
    os.makedirs(pathOut)

#load experiments
test_chr = Chromosome(name=chrin, centromere_search=False,
                      species='Homo sapiens', assembly='na')#, max_tad_size=260000)
test_chr.add_experiment(flag,exp_type='Hi-C', identifier='GM128Rao', resolution=res,
                       norm_data=matPath)



exp = test_chr.experiments[flag]
#this is for Hic PCHIC will be filter_julen and no norm
#exp.filter_columns(silent=False,draw_hist=False)
#exp.normalize_hic(silent=False, factor=None)
# If HiC data
if 'OneD' in matPath:
	exp.filter_columns(silent=False,draw_hist=False)
# If pcHiC or virtual pcHiC
else:
	PCHiC_filtering(exp)

#best parameters aFTER OPTIMISATION
#best cc
optpar = {'kforce': 5, 'lowfreq': low, 'upfreq': up, 'maxdist': maxd, 'scale': 0.01, 'reference': 'Merged'}

# I want the modelling to carsh if didnt stop in the time limit i said
#time1 = [int(t) for t in jobTime.split(':')]
#time1 = (time1[0] * 3600) + (time1[1] * 60) + time1[2] + (20 * 60) 
# will set time limit to the job limit in a way it crashes if not finished
if len(jobTime.split('-')) > 1:
    # fix jobtime
    hours = int(jobTime.split('-')[0]) * 24
    hours += int(jobTime.split('-')[1].split(':')[0])
    newtime = jobTime.split('-')[1].split(':')
    newtime[0] = str(hours)
    jobTime = ':'.join(newtime)
    time1 = [int(t) for t in jobTime.split(':')]
    time2 = (time1[0] * 3600) + (time1[1] * 60) + time1[2]
    time2 = time2 + (6 * 60)  # 5 min are the extra time before killing a job
else:
    time1 = [int(t) for t in jobTime.split(':')]
    time2 = (time1[0] * 3600) + (time1[1] * 60) + time1[2]
    time2 = time2 + (6 * 60)  # 5 min are the extra time before killing a job

# Build 3D models based on the HiC data. This is done by IMP.
#keep_restart_out_dir = path + 'lammpsSteps_mod/'
#if not os.path.exists(keep_restart_out_dir):
#    os.makedirs(keep_restart_out_dir)
#jobName = 'LF%sUF%sMdis%s_%sbp'%(str(low),str(up),str(maxd), str(res))
#keep_restart_out_dir = path + 'lammpsSteps_mod/jobArray_%s/' %jobName
#if not os.path.exists(keep_restart_out_dir):
#        os.makedirs(keep_restart_out_dir)
#restart_path = keep_restart_out_dir

keep_restart_out_dir = None
restart_path = False
models = exp.model_region(start=1,end=exp.size, n_models=nmodels, n_keep=nmodels,
    n_cpus=min(nmodels, n_cpu), 
    config=optpar,verbose=True, tool='lammps', 
	tmp_folder=lammpsOut, timeout_job=time2,
    cleanup=True, initial_conformation='random', 
    start_seed=random.choice(range(1000000)),
    hide_log=True,
    keep_restart_out_dir=keep_restart_out_dir,
    restart_path=restart_path,
    store_n_steps=2) #, connectivity='FENE')


save_models(models, pathOut+'%s_%s.modelsTemp'%(flag_name, lampsFlag))
print lammpsOut, pathOut+'%s_%s.modelsTemp'%(flag_name, lampsFlag)
#shutil.rmtree(keep_restart_out_dir) 
