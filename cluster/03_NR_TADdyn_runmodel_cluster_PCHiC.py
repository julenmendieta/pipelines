from pytadbit import Chromosome
from cPickle import dump, load
import sys
# CRG
sys.path.append('/users/mmarti/jmendieta/codigo')
sys.path.append('/users/mmarti/jmendieta/programas/PhD/3C')
# CNAG
sys.path.append('/home/devel/jmendietaesteban/PCHIC/codigo')
sys.path.append('/home/devel/jmendietaesteban/programas/PhD/3C')
import handy3C
import PCHiC
import argparse
import numpy as np
#import seaborn as sns
import matplotlib.pyplot as plt
import os, errno
import shutil
import random

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


parser = argparse.ArgumentParser(description='')
parser.add_argument('-l','--lowfreq', help='lowfreq',required=True)
parser.add_argument('-m','--maxdist',help='maxdist_range', required=True)
parser.add_argument('-d','--dcutoff_range',help='dcutoff_range', required=True)
parser.add_argument('-u','--upperfreq',help='upperfreq', required=True)
parser.add_argument('-lf','--lammpsfolder',help='folder_temp_lammps', required=True)
parser.add_argument('-p','--pathtomtrx',help='path_to_matrix', required=True)
parser.add_argument('-t','--jobtime',help='jobtime_HH:MM:SS', required=True)
parser.add_argument('-tp','--temp_path',help='path_to_tmp_files', required=False)

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

res = int(matPath.split('_')[-1][:-2])
lampsFlag = lammpsOut.split('_')[-1]
chrin='na'
path='/'.join(matPath.split('/')[:-1]) + '/'
res = int(matPath.split('_')[-1][:-2])
flag = matPath.split('/')[-2]

flag_name='%s_%sScaled01C%sL%sU%sM%sRes%s' %(flag, res, c, low, up, maxd, res)

# if we provide an alternative output dir we change it 
if tempOut != None:
    lammpsOut = tempOut


# Move to lammpsOut directory so its output gets stored there
#mkdir(lammpsOut)
#os.chdir(lammpsOut)

pathOut= path + 'finalModel/'
#matPath='/scratch/devel/jmendieta/Ferrer/modelling/PIdyn/reg32Norm/MatrixNormMin_43535000_44430000_compPI'

# Create outdir if does not exist
if not os.path.exists(pathOut):
    os.makedirs(pathOut)

#load experiments
test_chr = Chromosome(name=chrin, centromere_search=False,
                      species='Homo sapiens', assembly='na')#, max_tad_size=260000)
# If interaction index and not matrix
if matPath.split('/')[-1].startswith('interaction'):
        regionStart = int(matPath.split('/')[-1].split('_')[3].split('-')[1])
        regionEnd = int(matPath.split('/')[-1].split('_')[3].split('-')[2])
        matPath = handy3C.tableToMatrix(matPath, regionStart, regionEnd, resol=res, dirout="")
        test_chr.add_experiment(flag,exp_type='Hi-C', resolution=res,
                        norm_data=[matPath])
else:
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
	PCHiC.PCHiC_filtering(exp)

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
models = exp.model_region(start=1,end=exp.size, n_models=1, n_keep=1, n_cpus=1, config=optpar,verbose=True, tool='lammps', 
	tmp_folder=lammpsOut, timeout_job=time2,
    cleanup=True, initial_conformation='random', 
    start_seed=random.choice(range(1000000)),
    hide_log=True,
    keep_restart_out_dir=keep_restart_out_dir,
    restart_path=restart_path,
    store_n_steps=2) #, connectivity='FENE')


models.save_models(pathOut+'%s_%s.models'%(flag_name, lampsFlag))
print lammpsOut, pathOut+'%s_%s.models'%(flag_name, lampsFlag)
#shutil.rmtree(keep_restart_out_dir) 
