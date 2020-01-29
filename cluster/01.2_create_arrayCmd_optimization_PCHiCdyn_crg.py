import numpy as np
import time
import argparse

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

######################  Work with input  #####################
#lowfreq_arange = np.arange(-1.5,0.5,0.5)
#c_range=np.arange(200,500,100) #Cutoff
# El c_range optimo se saca asi
# resol * 2 * scale
#m_range=np.arange(200,500 ,100)
#upfreq_range=np.arange(-1, 1.5, 0.5)

parser = argparse.ArgumentParser(description='')
parser.add_argument('-lr','--lowfreqRange', help='lowfreqRange',required=True)
parser.add_argument('-cr','--cutoffRange', help='cutoffRange',required=True)
parser.add_argument('-mr','--maxdistRange',help='maxdist_range', required=True)
parser.add_argument('-ur','--upfreqRange', help='upfreqRange',required=True)
parser.add_argument('-p','--pathtomtrx',help='path_to_matrix', required=True)
parser.add_argument('-t','--jobtime',help='jobtime_HH:MM:SS', required=True)
parser.add_argument('-s','--scriptspath',help='path_to_scripts', required=True)
parser.add_argument('-nm','--nmodels',help='nmodels', required=True)

args = parser.parse_args()
lowfreq= [float(s) for s in cToDot(args.lowfreqRange).split('_')]
uperfreq= [float(s) for s in cToDot(args.upfreqRange).split('_')]
maxdist= [float(s) for s in cToDot(args.maxdistRange).split('_')]
dcutoff_range= args.cutoffRange
matPath=args.pathtomtrx
jobTime=args.jobtime
scriptsPath=args.scriptspath
nmodels=int(args.nmodels)

lowfreq_arange = np.arange(lowfreq[0],lowfreq[1],lowfreq[2])

# get dcutoff range to test probelms with maxdist
c_range = [float(s) for s in cToDot(args.cutoffRange).split('_')]


# El c_range optimo se saca asi
# resol * 2 * scale
m_range=np.arange(maxdist[0],maxdist[1] ,maxdist[2])
upfreq_range=np.arange(uperfreq[0],uperfreq[1], uperfreq[2])

flag = '_'.join(matPath.split('/')[-1].split('_')[1:3])
path='/'.join(matPath.split('/')[:-1]) + '/'
#jobTime = '04:00:00'  # HH:MM:SS
#path = "/scratch/devel/jmendieta/deLaat/modelling/reg4_WPL-KOD/"

#####################################################

print upfreq_range
print lowfreq_arange
print dcutoff_range
print m_range
# Keep record of the number of jobs to run
#njobs = len(lowfreq_arange)*len(c_range)*len(m_range)*len(upfreq_range)
#print njobs

fecha = time.strftime("%d-%m-%Y")
# File were we are going to keep all the posibilities of LF, UF etc to be run
runfile = '%s%s_%s.array'%(path, fecha, flag)
# Create a counter for the different lammps output directories
countl = 1
## Create the file with the commands to be run in the array
fout=open(runfile,'w')
for x in lowfreq_arange:
    #for c in c_range:
    for m in m_range:
        # check if we have a dcutoff small enough compared to maxdist to allow running
        # we allow an overlap of 50nm between maxdist and dcutoff
        if m - min(c_range) >= -50:
            for u in upfreq_range:
                if u >= x:
                    cmd= ''
                    cmd+= '%s01.2_NR_TADdyn_test_arg_PCHiC.py -l %s '%(scriptsPath, x)
                    cmd+= '-d %s '%dcutoff_range
                    cmd+= '-m %s '%m
                    cmd+= '-u %s '%u
                    cmd+= '-p %s '%matPath
                    cmd+= '-t %s '%jobTime
                    cmd+= '-nm %s '%str(nmodels)
                    cmd+= '\n'
                    fout.write(cmd)
                    countl += 1
fout.close()
njobs = countl - 1
print njobs, 'Jobs will be created'
## Create the file to launch the array

# asses priority in case of CRG cluster
if len(jobTime.split('-')) > 1:
    prior = 'long-sl7,mem_256,mem_512'
else:
    if int(jobTime.split(':')[0]) <= 6:
        prior = 'short-sl7,long-sl7,mem_256,mem_512,guest'
    else:
        prior = 'long-sl7,mem_256,mem_512'

fout=open('%s/arrayjobs.cmd'%(path),'w')
cmd=''
cmd+='''#!/bin/bash

#$ -N %s
#$ -o %s/$JOB_NAME_$JOB_ID_$TASK_ID.out
#$ -e %s/$JOB_NAME_$JOB_ID_$TASK_ID.err
#$ -t 1-%s
#$ -l h_rt=%s,virtual_free=16G
#$ -q %s
#$ -pe smp %s

module purge
module load GCCcore/6.3.0

# File were we have located our array commands
file=%s

# Get each command from the file and run them with python
orden=`sed "${SGE_TASK_ID}q;d" $file`
# will add the command for the temporal folder
orden=`echo $orden -tp $TMPDIR`

python $orden''' %(flag, path, path, njobs, jobTime, prior,
str(min(nmodels, 8)), runfile)

fout.write(cmd)
