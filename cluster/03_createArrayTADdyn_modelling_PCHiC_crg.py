import numpy as np
import time
import argparse

######################  TO CHANGE #####################
#nmodels = 500


#flag="Mod32"
#path = "/scratch/devel/jmendieta/Ferrer/modelling/PIdyn/reg32Norm/scale01"
#jobTime = '02:00:00'  # HH:MM:SS

#lowfreq=0
#uperfreq=0
#maxdist=400
#dcutoff=300
#res=5000
#lammpsOut="/scratch_tmp/%s/temp" %(flag)

#####################################################
parser = argparse.ArgumentParser(description='')
parser.add_argument('-l','--lowfreq', help='lowfreq',required=True)
parser.add_argument('-m','--maxdist',help='maxdist_range', required=True)
parser.add_argument('-d','--dcutoff_range',help='dcutoff_range', required=True)
parser.add_argument('-u','--upperfreq',help='upperfreq', required=True)
#parser.add_argument('-lf','--lammpsfolder',help='folder_temp_lammps', required=True)
parser.add_argument('-p','--pathtomtrx',help='path_to_matrix', required=True)
parser.add_argument('-t','--jobtime',help='jobtime_HH:MM:SS', required=True)
parser.add_argument('-s','--scriptspath',help='path_to_scripts', required=True)
parser.add_argument('-nm','--nmodels',help='output_nmodels', required=True)

args = parser.parse_args()
lowfreq=float(args.lowfreq)
uperfreq=float(args.upperfreq)
maxdist=float(args.maxdist)
dcutoff=float(args.dcutoff_range)
#lammpsOut=args.lammpsfolder
matPath=args.pathtomtrx
jobTime=args.jobtime
scriptsPath=args.scriptspath
nmodels=int(args.nmodels)


flag = matPath.split('/')[-2]
path='/'.join(matPath.split('/')[:-1]) + '/'
lammpsOut="/scratch_tmp/%s/temp" %(flag)



# Keep record of the number of jobs to run
njobs = nmodels
print njobs

fecha = time.strftime("%d-%m-%Y")
# File were we are going to keep all the posibilities of LF, UF etc to be run
runfile = '%s%s_%s.Modsarray'%(path, fecha, flag)
## Create the file with the commands to be run in the array
fout=open(runfile,'w')
for x in range(1, njobs + 1):
	cmd=''
	cmd+='%s03_NR_TADdyn_runmodel_cluster_PCHiC.py -l %s '%(scriptsPath, lowfreq)
	cmd+= '-d %s '%dcutoff
	cmd+= '-m %s '%maxdist
	cmd+= '-u %s '%uperfreq
	cmd+= '-lf %s '%(lammpsOut + '_' + str(x))
	cmd+= '-p %s '%matPath
	#cmd+= '-t %s '%jobTime # will add it in the arraMod.jobs so we can update it later
	cmd+= '\n'
	fout.write(cmd)
fout.close()


# asses priority in case of CRG cluster
if len(jobTime.split('-')) > 1:
        prior = 'long-sl7,mem_256,mem_512'
else:
        if int(jobTime.split(':')[0]) <= 6:
                prior = 'short-sl7,long-sl7,mem_256,mem_512,guest'

## Create the file to launch the array (sbatch file)
fout=open('%s/arrayModjobs.cmd'%(path),'w')
cmd=''
cmd+='''#!/bin/bash

#$ -N %s
#$ -o %s/$JOB_NAME_$JOB_ID_$TASK_ID.out
#$ -e %s/$JOB_NAME_$JOB_ID_$TASK_ID.err
#$ -t 1-%s
#$ -l h_rt=%s,virtual_free=16G
#$ -q %s
#$ -pe smp 1

module purge
module load GCC/5.3.0

# File were we have located our array commands
file=%s

# Get each command from the file and run them with python
orden=`sed "${SGE_TASK_ID}q;d" $file`
# will add the command for the temporal folder
orden=`echo $orden -tp $TMPDIR -t %s -jid ${SGE_TASK_ID}`

python $orden''' %(flag, path, path, njobs, jobTime, prior, runfile,
        jobTime)

fout.write(cmd)
