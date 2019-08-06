import os
import subprocess
import time
import copy
import numpy as np


# Function to turn the dots from a text into a c
# Useful to pass floating values into parser
def dotToC(text):
    newText = ''
    for t in text:
        if t == '.':
            newText += 'c'
        elif t == '-':
            newText += 'n'
        else:
            newText += t
    return newText

def optimPlots(scriptsPath, matPath):
        path='/'.join(matPath.split('/')[:-1]) + '/'
        # Join txt files
        script = scriptsPath + '02.1a_optTXTparser.py'
        cmd = script + ' -p %s -j yes' %path
        os.system("python %s" %cmd)

        # Create optimisation Plots
        script = scriptsPath + '02.1b_optimizationPlot_PCHiC.py'
        cmd = script + ' -p %s' %matPath
        os.system("python %s" %cmd)

def waitForJobsToFinish(clusterUser, maxJobs, cluster='cnag'):
    # check in which cluser we are
    if cluster == 'crg':
        clusterCheck = 'qstat'
        separateArrays='-g d'
        hideHeader=''
    elif cluster == 'cnag':
        clusterCheck = 'mnq'
        separateArrays='--array'
        hideHeader='--noheader'

    continueChecking = True
    while continueChecking == True:
        if cluster == 'crg':
            p = subprocess.Popen([clusterCheck, '%s -u %s' %(separateArrays, clusterUser)],
                                        shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = p.communicate()
            # in this case remove header from list
            jobs = list(line.split()[0] for line in out.split('\n')[2:-1])
        else:
            p = subprocess.Popen([clusterCheck, separateArrays, '-u %s' %clusterUser,
                                    hideHeader], 
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = p.communicate()
            jobs = list(line.split()[0] for line in out.split('\n')[:-1])
        # get all runing and about to be ran jobs
        runningJobs = len(jobs)
        if runningJobs > maxJobs:
            # Make some time 5 mins
            time.sleep(300)
        else:
            continueChecking = False

def saveDoneJobs(clusterUser, outTrash, jobPathsOptim, cluster='cnag'):
    if cluster == 'crg':
        clusterCheck = 'qstat'
    elif cluster == 'cnag':
        clusterCheck = 'mnq'

    # First get unfinished jobs
    if cluster == 'crg':
        p = subprocess.Popen([clusterCheck, '-u %s' %(clusterUser)],
                                shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
        # in this case remove header from list
        jobs = set(int(line.split()[0]) for line in out.split('\n')[2:-1])
    elif cluster == 'cnag':
        p = subprocess.Popen([clusterCheck, '-u %s' %clusterUser], stdout=subprocess.PIPE, 
                        stderr=subprocess.PIPE)
        out, err = p.communicate()
        jobs = set(int(line.split()[0].split('_')[0]) for line in out.split('\n') 
                        if line != '' and line.split()[0] != 'JOBID')
    # Second we load prevoius file if exists
    toAdd = []
    if os.path.exists(outTrash + 'optimList.txt'):
        with open(outTrash + 'optimList.txt', 'r') as f:
            for line in f:
                if len(line) != 0:
                    line = line.split()
                    line = [line[0], int(line[1])]
                    # check if some jobs didnt finish
                    if int(line[1]) in jobs:
                        toAdd.append(line)
    # Merge both lists
    jobPathsOptim += toAdd

    # Recreate the file
    with open(outTrash + 'optimList.txt', 'w') as f:
        for j in jobPathsOptim:
            f.write('%s\t%s\n' %(j[0], j[1]))

def createFailedOptimeRunFiles(matPath, combinations, scriptsPath, jobTime, nmodelsOptim,
                                prior='long-sl7,mem_256,mem_512', cluster=False):
    flag = matPath.split('/')[-2]
    path='/'.join(matPath.split('/')[:-1]) + '/'

    fecha = time.strftime("%d-%m-%Y")
    runfile = '%s%s_%s.array'%(path, fecha, flag)
    # Create a counter for the different lammps output directories
    countl = 1
    ## Create the file with the commands to be run in the array
    fout=open(runfile,'w')
    for c,m,l,u in combinations.keys():
        cmd=''
        cmd+='%s01.2_NR_TADdyn_test_arg_PCHiC.py -l %s '%(scriptsPath, l)
        cmd+= '-d %s '%c
        cmd+= '-m %s '%m
        cmd+= '-u %s '%u
        cmd+= '-p %s '%matPath
        cmd+= '-t %s '%jobTime
        cmd+= '-nm %s '%str(nmodelsOptim)
        cmd+= '\n'
        fout.write(cmd)
        countl += 1
    fout.close()
    njobs = countl - 1
    print njobs, 'Jobs will be created'
    ## Create the file to launch the array
    fout=open('%s/arrayjobs.cmd'%(path),'w')
    cmd=''
    # different cases when CRG or CNAG cluster
    if cluster == 'crg':
        cmd+='''#!/bin/bash

#$ -N %s
#$ -o %s/$JOB_NAME_$JOB_ID_$TASK_ID.out
#$ -e %s/$JOB_NAME_$JOB_ID_$TASK_ID.err
#$ -t 1-%s
#$ -l h_rt=%s,virtual_free=16G
#$ -q %s
#$ -pe smp %s

module purge
module load GCC/5.3.0

# File were we have located our array commands
file=%s

# Get each command from the file and run them with python
orden=`sed "${SGE_TASK_ID}q;d" $file`
# will add the command for the temporal folder
orden=`echo $orden -tp $TMPDIR`
echo $orden
python $orden''' %(flag, path, path, njobs, jobTime, prior, str(min(nmodelsOptim, 8)), runfile)

    elif cluster == 'cnag':
        cmd+='''#!/bin/bash

#SBATCH --job-name=%s
#SBATCH --output=%s/%%A_submatrix_%%a.out
#SBATCH --error=%s/%%A_submatrix_%%a.err
#SBATCH --array=1-%s%%300
#SBATCH --time=%s
#SBATCH --qos=%s
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=%s

module purge

# File were we have located our array commands
file=%s

# Get each command from the file and run them with python
orden=`sed "${SLURM_ARRAY_TASK_ID}q;d" $file`
python $orden''' %(flag, path, path, njobs, jobTime, prior, str(min(nmodelsOptim, 8)), runfile)
    else:
        No_contemplado

    fout.write(cmd)


def cleaner(matPath, justClean=False):
    path = '/'.join(matPath.split('/')[:-1]) + '/'
    toCompress = False

    allfiles = os.listdir(path)
    # Remove all job status files
    for a in allfiles:
        if a.endswith('err'):
            p = subprocess.Popen(['rm', '%s%s' %(path, a)],
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = p.communicate()
        elif a.endswith('out'):
            p = subprocess.Popen(['rm', '%s%s' %(path, a)],
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = p.communicate()

        # Check presence of opt_ files
        if a.startswith('opt_LF'):
            toCompress = True

    # compress all optimization output and model files
    if toCompress == True and justClean == False:
        _ = os.system('tar -cvf %soptimOut.tar %sopt_LF* --remove-files' %(path, path))



# Select what to prepare
jobPathsOptim = []
jobPathsModel = []
# Get date from when we ran this script
date = time.strftime("%Y-%m-%d")


####################### TO CHANGe ###################s####################
# create optimisation files
optimization = False
# run optimisation
runOptim = False
# join optimisation files and make plots
optimOut = False
# rerun the ones which failed (checks from matPaths)
runOptimFailed = False
cleanFoldersOptim = False
if runOptimFailed:
    # dcutoff list if you want just to rerun the ones from here
    # To be writen as range (beggin, end, step)
    c_focus = [250.0,450.0,50.0]
    # maxdist list if you want just to rerun the ones from here
    # To be writen as range (beggin, end, step)
    m_focus = [200.0,700.0, 100.0]
# this will create a modelling file with the parameters combinations with top correlation
# MUST be checked for cases of combinations like low -1 and up 1
createModellingFile = False
if createModellingFile:
    # to check distirbution of correlations at each dcutoff
    show_dcut = False  # False to not use, True to check that dcutoff
    # to check distribution of correlations at certain dcutoff separating by maxdist
    dcut = False  # False not to check, integer or float with dcutoff to do it
    # create file with maximum correlations at certain dcutoff and maxdist
    # dcutof and maxdist parameters where to get the top correlators
    dcut_ = 200
    maxd_ = 400
    jobTime2_ = '03:00:00'

modelling = True
joinModels = False
checkModellinTime = False
cleanFolders= False

jobTime = '0-03:00:00'
# CNAG
scriptsPath = '/scratch/devel/jmendieta/PCHiC_assesment/modelling/'
outTrash = '/scratch/devel/jmendieta/PCHiC_assesment/modelling/'
loadMatrixFile = '/scratch/devel/jmendieta/PCHiC_assesment/modelling/matrixList.txt'
cluster = 'cnag'
# CRG
#scriptsPath = '/users/mmarti/jmendieta/scratch/PCHiC_assesment/modelling/'
#outTrash = '/users/mmarti/jmendieta/scratch/PCHiC_assesment/modelling/'
#loadMatrixFile = '/users/mmarti/jmendieta/scratch/PCHiC_assesment/modelling/matrixList.txt'
#cluster = 'crg'

# You name in the cluster user queue
clusterUser = 'jmendieta'
# Maximum number of jobs allowed to be running before launching another job array
maxJobs = 1000
# Get matrices from file
matPaths= []
with open(loadMatrixFile,'r') as f:
    for line in f:
        matPaths += line.split()

nmodelsOptim = 100
lowfreq_arange = [-1.0,0.5,0.5]
c_range= [50,450,50] #Cutoff no more overlap than 50nm with maxdist will be allowed
# El c_range optimo es menor o igual que la m_range optima,
# que se saca asi
# resol *scale * 2
# (resol * scale = radius) * 2 = diameter of a particle
m_range= [100,500,100]
upfreq_range= [-1,1.0,0.5]


####################### END CHANGE ####################################


if cluster == 'crg':
    runArraycmd = 'qsub'
    clusterCheck = 'qstat'
    separateArrays='-g d'
    hideHeader=''
    # asses priority
    if len(jobTime.split('-')) > 1:
        prior = 'long-sl7,mem_256,mem_512'
        # fix jobtime
        hours = int(jobTime.split('-')[0]) * 24
        hours += int(jobTime.split('-')[1].split(':')[0])
        newtime = jobTime.split('-')[1].split(':')
        newtime[0] = str(hours)
        jobTime = ':'.join(newtime)
    else:
        if int(jobTime.split(':')[0]) <= 6:
            prior = 'short-sl7,long-sl7,mem_256,mem_512,guest'
elif cluster == 'cnag':
    runArraycmd = 'sbatch'
    clusterCheck = 'mnq'
    separateArrays='--array'
    hideHeader='--noheader'
    prior = 'lowprio'
#matPath = '/scratch/devel/jmendieta/deLaat/modelling/reg4_WPL-KOD/MatrixFreqNorm_WPL-KOD_reg4_chr8-120780000-122030000_10000bp'
#matPaths = ['/scratch/devel/jmendieta/PCHiC_assesment/modelling/PI/reg1/MatrixNormFreqFiltrd_PI_reg1_chr8-117680000-118410000_5000bp']

## Define optimization parameters
if optimization == True:
    for matPath in matPaths:
        script = scriptsPath + '01.2_create_arrayCmd_optimization_PCHiCdyn_%s.py' %(cluster)
        lowfreq_arangeDot = dotToC('_'.join(str(j) for j in lowfreq_arange)) 
        c_rangeDot= dotToC('_'.join(str(j) for j in c_range))
        m_rangeDot= dotToC('_'.join(str(j) for j in m_range))
        upfreq_rangeDot= dotToC('_'.join(str(j) for j in upfreq_range))

        cmd = script + ' -lr %s -cr %s -mr %s -ur %s -p %s -t %s -s %s -nm %s' %(lowfreq_arangeDot,
                                c_rangeDot,
                        m_rangeDot,
                        upfreq_rangeDot,
                        matPath,
                        jobTime,
                        scriptsPath,
                        str(nmodelsOptim))
        print cmd
        os.system("python %s" %cmd)

        if runOptim == True:
            # We will just allow $maxJobs jobs to be running more or less
            waitForJobsToFinish(clusterUser, maxJobs, cluster=cluster)
            
            # run batch job and keep job id
            path='/'.join(matPath.split('/')[:-1]) + '/'
            p = subprocess.Popen([runArraycmd, '%sarrayjobs.cmd' %path], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = p.communicate()
            if cluster == 'crg':
                jobId = out.split()[2].split('.')[0]
            elif cluster == 'cnag':
                jobId = out.split()[-1]    
            print jobId
            # assign job to current path
            jobPathsOptim.append([matPath, int(jobId)])

        ## Save all jobs just in case
        saveDoneJobs(clusterUser, outTrash, jobPathsOptim, cluster=cluster)
    


if optimOut == True:
    # First check if we have optimisation jobs to prepare output
    if os.path.exists(outTrash + 'optimList.txt'):
        with open(outTrash + 'optimList.txt', 'r') as f:
            for line in f:
                line = line.split()
                line = [line[0], int(line[1])]
                # check if some jobs didnt finish
                if line not in jobPathsOptim:
                    jobPathsOptim.append(line)


    # Keep record of the number of jobs we finished
    doneJobs = 0
    jobs = range(3)
    # While we still have jobs
    finish = False
    jobPathsOptim_ = copy.deepcopy(jobPathsOptim)
    while finish == False:
        # Check jobs that are still running
        if cluster == 'crg':
            p = subprocess.Popen([clusterCheck, '%s -u %s' %(separateArrays, clusterUser)],
                        shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = p.communicate()
            # in this case remove header from list
            jobs = set(int(line.split()[0]) for line in out.split('\n')[2:-1])
        elif cluster == 'cnag':
            p = subprocess.Popen([clusterCheck, separateArrays, '-u %s' %clusterUser, 
                hideHeader], 
                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = p.communicate()
            jobs = set(int(line.split()[0].split('_')[0]) for line in out.split('\n')[:-1])
        # check which of the jobs from the optimisation are still running
        oneFound = False
        toKeep = []
        for nk, kjob in enumerate(jobPathsOptim_):
            jobId = kjob[1]
            matPath = kjob[0]
            if jobId not in jobs:
                optimPlots(scriptsPath, matPath)
            else:
                toKeep.append(nk)
                oneFound = True
        jobPathsOptim_ = [jobPathsOptim_[k] for k in toKeep]
        # If we finished all the jobs from the list but there are others around
        if oneFound == False:
            finish = True 

    # To finish we empty the optimisation josbs file
    with open(outTrash + 'optimList.txt', 'w') as f:
        pass

## If want to rerun the ones which failed
if runOptimFailed == True:
    ###################### TO CHANGE #######################
    # dcutoff list if you want just to rerun the ones from here
    # To be writen as range (beggin, end, step)
    #c_focus = [200.0, 300.0, 100.0]
    # maxdist list if you want just to rerun the ones from here
    # To be writen as range (beggin, end, step)
    #m_focus = [300.0,400.0, 100.0]
    ###################### END-CHANGE #######################

    for matPath in matPaths:
        fi0 = '/'.join(matPath.split('/')[:-1]) + '/'

        # Check if we want to focus somewhere and select range accordingly
        if c_focus != False:
            cRange = np.arange(c_focus[0], c_focus[1], c_focus[2])
            cRtext = dotToC('_'.join(str(j) for j in c_focus))
        else:
            cRange = np.arange(c_range[0], c_range[1], c_range[2])
            cRtext = dotToC('_'.join(str(j) for j in c_range))
        if m_focus != False:
            mRange = np.arange(m_focus[0], m_focus[1], m_focus[2])
        else:
            mRange = np.arange(m_range[0], m_range[1], m_range[2])

        # get all combinations that should be
        combinations = {}
        for c in cRange:
            for m in mRange:
                # we allow an overlap of 50nm
                if (m - c) >= -50:
                    for l in np.arange(lowfreq_arange[0], lowfreq_arange[1], lowfreq_arange[2]):
                        for u in np.arange(upfreq_range[0], upfreq_range[1], upfreq_range[2]):
                            if u >= l:
                                combinations[(c,m,l,u)] = False

        # Mark as True the ones we actually have
        for c in cRange:
            fi = fi0 + 'ValOptimisationC%s.txt' %float(c)
            # If at least one job entered
            if os.path.exists(fi):
                with open(fi, 'r') as f:
                    for line in f:
                        if not line.startswith('#'):
                            line = line.split()
                            combi = (int(line[5]), int(line[2]), float(line[3]), float(line[4]))
                            combinations[combi] = True
            else:
                print 'WARNING!!!'
                print 'None of the jobs finished in:\n%s' %fi
                print 'This is very unlikely, re-running, but check everything is ok'
                print 'WARNING!!!'

        # Keep just the ones we didnt got
        for key, val in combinations.items():
            if val == True:
                del combinations[key]

        # If we have areas to re-do
        if len(combinations.keys()):
            # re-run failed jobs
            # We will just allow $maxJobs jobs to be running more or less
            waitForJobsToFinish(clusterUser, maxJobs, cluster=cluster)

            # merge all dcutoff variations
            combinations2 = {}
            for comb in combinations.keys():
                combinations2[tuple([cRtext] + list(comb[1:]))] = True


            # Modify job files
            path='/'.join(matPath.split('/')[:-1]) + '/'
            createFailedOptimeRunFiles(matPath, combinations2, scriptsPath, jobTime, nmodelsOptim,
                    prior=prior, cluster=cluster)

            # run batch job and keep job id
            p = subprocess.Popen([runArraycmd, '%sarrayjobs.cmd' %path], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = p.communicate()
            if cluster == 'crg':
                jobId = out.split()[2].split('.')[0]
            elif cluster == 'cnag':
                jobId = out.split()[-1]
            print jobId
            # assign job to current path
            jobPathsOptim.append([matPath, int(jobId)])

    ## Save all jobs just in case
    saveDoneJobs(clusterUser, outTrash, jobPathsOptim, cluster=cluster)



if modelling == True:
    # here each job is one model, before each job was 100 in chunks of 8
    maxJobs = maxJobs * 8

    # Check if we filled the file to set the modelling parameters
    toModel = []
    if os.path.exists(outTrash + 'modellinParams.txt'):
        # Open model characteristics file
        with open(outTrash + 'modellinParams.txt', 'r') as f:
            for line in f:
                line = line.split()
                # check if some jobs didnt finish
                toModel.append(line)
        
        # get print when same matrix path (array file will be overwriten)
        oldPath1 = ''
        for tm in toModel:
            # if empty line
            if len(tm) == 0:
                continue
            # toModel seria un fichero con
            # matPath lowfreq upfreq d_cutoff maxdist jobTime nmodels
            matPath = tm[0]
            lowfreq = float(tm[1])
            upfreq = float(tm[2])
            d_cutoff= float(tm[3])
            # El c_range optimo es menor o igual que la m_range optima,
            # que se saca asi
            # resol * 2 * scale
            maxdist= float(tm[4])
            jobTime = tm[5]  # HH:MM:SS
            nmodels = int(tm[6])

            # warn if same matPath as previous
            if oldPath1 == matPath:
                print 'you are modelling the same matrix in same folder, arrayjobs\
 file will be overwriten, so just last job in folder will be done many times'
            oldPath1 = matPath
            ##################### TO CHANGE ############################################
            #lowfreq = 0.0
            #d_cutoff= 300 #Cutoff
            # El c_range optimo es menor o igual que la m_range optima,
            # que se saca asi
            # resol * 2 * scale
            #maxdist= 400
            #upfreq = 0.5
            #jobTime = '02:00:00'  # HH:MM:SS
            #nmodels = 500
            ####################### END CHANGE ####################################
            script = scriptsPath + '03_createArrayTADdyn_modelling_PCHiC_%s.py' %(cluster)

            cmd = script + ' -l %s -d %s -m %s -u %s -p %s -t %s -s %s -nm %s' %(lowfreq,
                                            d_cutoff,
                                        maxdist,
                                        upfreq,
                                        matPath,
                                        jobTime,
                                        scriptsPath,
                                nmodels)
            print cmd
            os.system("python %s" %cmd)

            # We will just allow $maxJobs jobs to be running more or less
            waitForJobsToFinish(clusterUser, maxJobs, cluster=cluster)

            # run batch job and keep job id
            path='/'.join(matPath.split('/')[:-1]) + '/'
            p = subprocess.Popen([runArraycmd, '%sarrayModjobs.cmd' %path], 
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = p.communicate()
            if cluster == 'crg':
                jobId = out.split()[2].split('.')[0]
            elif cluster == 'cnag':
                jobId = out.split()[-1]
            print jobId
            # assign job to current path
            jobPathsModel.append([matPath, jobId])
            print jobPathsModel

        # Store all launched jobs
        # First get unfinished jobs
        if cluster == 'crg':
            p = subprocess.Popen([clusterCheck, '-u %s' %(clusterUser)],
                    shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = p.communicate()
            # in this case remove header from list
            jobs = set(int(line.split()[0]) for line in out.split('\n')[2:-1])
        else:
            p = subprocess.Popen([clusterCheck, '-u %s' %clusterUser], stdout=subprocess.PIPE, 
                    stderr=subprocess.PIPE)
            out, err = p.communicate()
            jobs = set(int(line.split()[0].split('_')[0]) for line in out.split('\n') 
                    if line != '' and line.split()[0] != 'JOBID')
        # Second we load prevoius file if exists
        toAdd = []
        if os.path.exists(outTrash + 'modelList.txt'):
            with open(outTrash + 'modelList.txt', 'r') as f:
                for line in f:
                    line = line.split()
                    if len(line) != 0:
                        # check if some jobs didnt finish
                        if int(line[1]) in jobs:
                            toAdd.append([line[0], int(line[1])])
        # Merge both lists
        jobPathsModel += toAdd

        # Recreate the file
        with open(outTrash + 'modelList.txt', 'w') as f:
            for j in jobPathsModel:
                f.write('%s\t%s\n' %(j[0], j[1]))

    else:
        print 'No modelling parameters file found in:\n%s' %(outTrash + 'modellinParams.txt')

## Join models
########### No esta para nada terminado
if joinModels == True:
    # First check if we have modelling jobs to prepare output
    toAdd = []
    if os.path.exists(outTrash + 'modelList.txt'):
        with open(outTrash + 'modelList.txt', 'r') as f:
            for line in f:
                line = line.split()
                line = [line[0], int(line[1])]
                # check if some jobs didnt finish
                if line not in jobPathsModel:
                    jobPathsModel.append(line)


    # Keep record of the number of jobs we finished
    doneJobs = 0
    # Get active jobs
    if cluster == 'crg':
        p = subprocess.Popen([clusterCheck, '-u %s' %(clusterUser)],
                    shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
        # in this case remove header from list
        jobs = set(int(line.split()[0]) for line in out.split('\n')[2:-1])
    else:
        p = subprocess.Popen([clusterCheck, '-u %s' %clusterUser], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
        jobs = set(int(line.split()[0].split('_')[0]) for line in out.split('\n') if line != '' and line.split()[0] != 'JOBID')

    
    # While we still have jobs
    finish = False
    while finish == False:
        # Check status of optimization
        for kjob in jobPathsModel:
            if cluster == 'crg':
                p = subprocess.Popen([clusterCheck, '-u %s' %(clusterUser)],
                        shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                out, err = p.communicate()
                # in this case remove header from list
                jobs = set(int(line.split()[0]) for line in out.split('\n')[2:-1])
            else:
                p = subprocess.Popen([clusterCheck, '-u %s' %clusterUser], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                out, err = p.communicate()
                jobs = set(int(line.split()[0].split('_')[0]) for line in out.split('\n') if line != '' and line.split()[0] != 'JOBID')
            # Mirar que siga solo cuando tengo jobs de la lista pendientes            
            # If finish
            jobId = kjob[1]
            if jobId not in jobs:
                doneJobs += 1
                matPath = kjob[0]
                # we merge the jobs
                mpath = '/'.join(matPath.split('/')[:-1]) + '/finalModel/'
                toRun = scriptsPath + '06_combineModels.py %s' %mpath    
                # Run script for merging
                print "python %s" %toRun
                os.system("python %s" %toRun)
                #optimPlots(scriptsPath, matPath)
        # If we finished all the jobs from the list but there are others around
        if doneJobs == len(jobPathsModel):
            finish = True 

    # To finish we empty the optimisation josbs file
    with open(outTrash + 'modelList.txt', 'w') as f:
        pass

# check if the time limit for modelling part is ok
# NOT OPTIMISED FOR CRG CLUSTER
if checkModellinTime == True:
    p = subprocess.Popen(['sacct', '-S %s -u %s --format=jobid,state,elapsed' %(date, clusterUser)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p = subprocess.check_output("sacct -S %s -u %s --format=jobid,state,elapsed" %(date, clusterUser), shell=True)


# clean files and group output
if cleanFolders == True:
    # First check if we have optimisation jobs to prepare output
    if os.path.exists(outTrash + 'modelList.txt'):
        with open(outTrash + 'modelList.txt', 'r') as f:
            for line in f:
                line = line.split()
                line = [line[0], int(line[1])]
                # check if some jobs didnt finish
                if line not in jobPathsOptim:
                    jobPathsOptim.append(line)


    # Keep record of the number of jobs we finished
    doneJobs = 0
    jobs = range(3)
    # While we still have jobs
    finish = False
    jobPathsOptim_ = copy.deepcopy(jobPathsOptim)
    while finish == False:
        # Check jobs that are still running
        if cluster == 'crg':
            p = subprocess.Popen([clusterCheck, '-u %s' %(clusterUser)],
                    shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = p.communicate()
            # in this case remove header from list
            jobs = set(int(line.split()[0]) for line in out.split('\n')[2:-1])
        else:
            p = subprocess.Popen([clusterCheck, '-u %s' %clusterUser], stdout=subprocess.PIPE, 
                    stderr=subprocess.PIPE)
            out, err = p.communicate()
            jobs = set(int(line.split()[0].split('_')[0]) for line in out.split('\n') 
                    if line != '' and line.split()[0] != 'JOBID')

        # check which of the jobs from the optimisation are still running
        oneFound = False
        toKeep = []
        for nk, kjob in enumerate(jobPathsOptim_):
            jobId = kjob[1]
            matPath = kjob[0]
            if jobId not in jobs:
                cleaner(matPath)
            else:
                toKeep.append(nk)
                oneFound = True
        jobPathsOptim_ = [jobPathsOptim_[k] for k in toKeep]
        # If we finished all the jobs from the list but there are others around
        if oneFound == False:
            finish = True

# Clean files after optimisation
if cleanFoldersOptim == True:
    # First check if we have optimisation jobs to prepare output
    if os.path.exists(outTrash + 'optimList.txt'):
        with open(outTrash + 'optimList.txt', 'r') as f:
            for line in f:
                line = line.split()
                line = [line[0], int(line[1])]
                # check if some jobs didnt finish
                if line not in jobPathsOptim:
                    jobPathsOptim.append(line)


    # Keep record of the number of jobs we finished
    doneJobs = 0
    jobs = range(3)
    # While we still have jobs
    finish = False
    jobPathsOptim_ = copy.deepcopy(jobPathsOptim)
    while finish == False:
        # Check jobs that are still running
        if cluster == 'crg':
            p = subprocess.Popen([clusterCheck, '-u %s' %(clusterUser)],
                    shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = p.communicate()
            # in this case remove header from list
            jobs = set(int(line.split()[0]) for line in out.split('\n')[2:-1])
        else:
            p = subprocess.Popen([clusterCheck, '-u %s' %clusterUser], stdout=subprocess.PIPE, 
                    stderr=subprocess.PIPE)
            out, err = p.communicate()
            jobs = set(int(line.split()[0].split('_')[0]) for line in out.split('\n') 
                    if line != '' and line.split()[0] != 'JOBID')

        # check which of the jobs from the optimisation are still running
        oneFound = False
        toKeep = []
        for nk, kjob in enumerate(jobPathsOptim_):
            jobId = kjob[1]
            matPath = kjob[0]
            if jobId not in jobs:
                cleaner(matPath, justClean=True)
            else:
                toKeep.append(nk)
                oneFound = True
        jobPathsOptim_ = [jobPathsOptim_[k] for k in toKeep]
        # If we finished all the jobs from the list but there are others around
        if oneFound == False:
            finish = True

    

    
if createModellingFile == True:
    topModels = '%s_%s_%s' %(dcut_, maxd_, jobTime2_)

    script = scriptsPath + '02.3_Check_optimizationOutput.py'
    cmd = script + ' -op %s -sdc %s -dc %s -tm %s' %(outTrash, str(show_dcut), str(dcut), topModels)
    os.system("python %s" %cmd)
