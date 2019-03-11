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
        script = scriptsPath + 'optTXTparser.py'
        cmd = script + ' -p %s' %path
        os.system("python %s" %cmd)

        # Create optimisation Plots
        script = scriptsPath + '02.1_optimizationPlot_PCHiC.py'
        cmd = script + ' -p %s' %matPath
        os.system("python %s" %cmd)

def waitForJobsToFinish(clusterUser, maxJobs):
	continueChecking = True
	while continueChecking == True:
		p = subprocess.Popen(['mnq', '--array', '-u %s' %clusterUser, '--noheader'], 
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

def saveDoneJobs(clusterUser, outTrash, jobPathsOptim):
	# First get unfinished jobs
	p = subprocess.Popen(['mnq', '-u %s' %clusterUser], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = p.communicate()
	jobs = set(int(line.split()[0].split('_')[0]) for line in out.split('\n') if line != '' and line.split()[0] != 'JOBID')
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

def createFailedOptimeRunFiles(matPath, combinations, scriptsPath, jobTime, nmodelsOptim):
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
    cmd+='''#!/bin/bash

#SBATCH --job-name=%s
#SBATCH --output=%s/%%A_submatrix_%%a.out
#SBATCH --error=%s/%%A_submatrix_%%a.err
#SBATCH --array=1-%s%%300
#SBATCH --time=%s
#SBATCH --qos=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=%s

module purge

# File were we have located our array commands
file=%s

# Get each command from the file and run them with python
orden=`sed "${SLURM_ARRAY_TASK_ID}q;d" $file`
python $orden''' %(flag, path, path, njobs, jobTime, str(min(nmodelsOptim, 8)), runfile)

    fout.write(cmd)

# Select what to prepare
jobPathsOptim = []
jobPathsModel = []
# Get date from when we ran this script
date = time.strftime("%Y-%m-%d")


################### Queda cambiar para que en RAOvirt cargue los modelos con la funcion para tablas
####################### TO CHANGe ###################s####################
optimization = True
runOptim = True
optimOut = False
runOptimFailed = False
modelling = False
joinModels = False
checkModellinTime = False

jobTime = '06:00:00'
scriptsPath = '/scratch/devel/jmendieta/PCHiC_assesment/modelling/'
outTrash = '/scratch/devel/jmendieta/PCHiC_assesment/modelling/'
loadMatrixFile = '/scratch/devel/jmendieta/PCHiC_assesment/modelling/matrixList.txt'
# You name in the cluster user queue
clusterUser = 'jmendieta'
# Maximum number of jobs allowed to be running before launching another job array
maxJobs = 800
# Get matrices from file
matPaths= []
with open(loadMatrixFile,'r') as f:
	for line in f:
		matPaths += line.split()

nmodelsOptim = 100
lowfreq_arange = [-1.5,0.5,0.5]
c_range= [200,300,100] #Cutoff
# El c_range optimo es menor o igual que la m_range optima,
# que se saca asi
# resol * 2 * scale
m_range= [200,500,100]
upfreq_range= [-1,1.5,0.5]

####################### END CHANGE ####################################
#matPath = '/scratch/devel/jmendieta/deLaat/modelling/reg4_WPL-KOD/MatrixFreqNorm_WPL-KOD_reg4_chr8-120780000-122030000_10000kb'
#matPaths = ['/scratch/devel/jmendieta/PCHiC_assesment/modelling/PI/reg1/MatrixNormFreqFiltrd_PI_reg1_chr8-117680000-118410000_5000bp']




## Define optimization parameters
if optimization == True:
	for matPath in matPaths:
		script = scriptsPath + '01.2_create_arrayCmd_optimization_PCHiCdyn.py'
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
			waitForJobsToFinish(clusterUser, maxJobs)
			
			# run batch job and keep job id
			path='/'.join(matPath.split('/')[:-1]) + '/'
			p = subprocess.Popen(['sbatch', '%sarrayjobs.cmd' %path], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			out, err = p.communicate()
			jobId = out.split()[-1]	
			print jobId
			# assign job to current path
			jobPathsOptim.append([matPath, int(jobId)])

	## Save all jobs just in case
	saveDoneJobs(clusterUser, outTrash, jobPathsOptim)
	


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
		p = subprocess.Popen(['mnq', '--array', '-u %s' %clusterUser, '--noheader'], 
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
	c_focus = [200, 300, 100]
	# maxdist list if you want just to rerun the ones from here
	# To be writen as range (beggin, end, step)
	m_focus = [300,500, 100]
	###################### END-CHANGE #######################

	for matPath in matPaths:
	    fi0 = '/'.join(matPath.split('/')[:-1]) + '/'

	    # Check if we want to focus somewhere and select range accordingly
	    if c_focus != False:
		cRange = np.arange(c_focus[0], c_focus[1], c_focus[2])
	    else:
		cRange = np.arange(c_range[0], c_range[1], c_range[2])

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
		waitForJobsToFinish(clusterUser, maxJobs)

		# Modify job files
		path='/'.join(matPath.split('/')[:-1]) + '/'
		createFailedOptimeRunFiles(matPath, combinations, scriptsPath, jobTime, nmodelsOptim)

		# run batch job and keep job id
		p = subprocess.Popen(['sbatch', '%sarrayjobs.cmd' %path], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		out, err = p.communicate()
		jobId = out.split()[-1]	
		print jobId
		# assign job to current path
		jobPathsOptim.append([matPath, int(jobId)])

	## Save all jobs just in case
	saveDoneJobs(clusterUser, outTrash, jobPathsOptim)



if modelling == True:
	lalalalal
	############### QuEDA PONER EL LIMITE DE JOBS QUE PUEDEN EJECUTARSE
	# Check if we filled the file to set the modelling parameters
	toModel = []
	if os.path.exists(outTrash + 'modellinParams.txt'):
		# Open model characteristics file
		with open(outTrash + 'modellinParams.txt', 'r') as f:
			for line in f:
				line = line.split()
				# check if some jobs didnt finish
				toModel.append(line)

		for tm in toModel:
			# toModel seria un fichero con
			# matPath lowfreq upfreq d_cutoff maxdist jobTime nmodels
			matPath = toModel[0]
			lowfreq = toModel[1]
			upfreq = toModel[2]
			d_cutoff= toModel[3]
			# El c_range optimo es menor o igual que la m_range optima,
			# que se saca asi
			# resol * 2 * scale
			maxdist= toModel[4]
			jobTime = toModel[5]  # HH:MM:SS
			nmodels = toModel[6]
			##################### TO CHANGE ############################################
			lowfreq = 0.0
			d_cutoff= 300 #Cutoff
			# El c_range optimo es menor o igual que la m_range optima,
			# que se saca asi
			# resol * 2 * scale
			maxdist= 400
			upfreq = 0.5
			jobTime = '02:00:00'  # HH:MM:SS
			nmodels = 500
			####################### END CHANGE ####################################
			script = scriptsPath + '03_createArrayTADdyn_modelling_PCHiC.py'

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

			# run batch job and keep job id
			path='/'.join(matPath.split('/')[:-1]) + '/'
			p = subprocess.Popen(['sbatch', '%sarrayModjobs.cmd' %path], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			out, err = p.communicate()
			jobId = int(out.split()[-1])
			print jobId
			# assign job to current path
			jobPathsModel.append([matPath, jobId])
			print jobPathsModel

		# Store all launched jobs
		# First get unfinished jobs
		p = subprocess.Popen(['mnq', '-u jmendieta'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		out, err = p.communicate()
		jobs = set(int(line.split()[0].split('_')[0]) for line in out.split('\n') if line != '' and line.split()[0] != 'JOBID')
		# Second we load prevoius file if exists
		toAdd = []
		if os.path.exists(outTrash + 'modelList.txt'):
			with open(outTrash + 'modelList.txt', 'r') as f:
				for line in f:
					line = line.split()
					# check if some jobs didnt finish
					if int(line[1]) in jobs:
						toAdd.append([line[0], int(line[1])])
		# Merge both lists
		jobPathsModel += toAdd

		# Recreate the file
		with open(outTrash + 'modelList.txt', 'w') as f:
			for j in jobPathsOptim:
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
					jobPathsModel += line


	# Keep record of the number of jobs we finished
	doneJobs = 0
	# Get active jobs
	p = subprocess.Popen(['mnq', '-u jmendieta'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = p.communicate()
	jobs = set(int(line.split()[0].split('_')[0]) for line in out.split('\n') if line != '' and line.split()[0] != 'JOBID')

	
	# While we still have jobs
	finish = False
	while len(jobs) > 0 and finish == False:
		# Check status of optimization
		for kjob in jobPathsModel:
			p = subprocess.Popen(['mnq', '-u jmendieta'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			out, err = p.communicate()
			# Mirar que siga solo cuando tengo jobs de la lista pendientes			
			# If finish
			jobs = set(int(line.split()[0].split('_')[0]) for line in out.split('\n') if line != '' and line.split()[0] != 'JOBID')
			jobId = kjob[1]
			if jobId not in jobs:
				doneJobs += 1
				matPath = kjob[0]
				optimPlots(scriptsPath, matPath)
		# If we finished all the jobs from the list but there are others around
		if doneJobs == len(jobPathsModel):
			finish = True 



	# To finish we empty the optimisation josbs file
	with open(outTrash + 'modelList.txt', 'w') as f:
		pass

# check if the time limit for modelling part is ok
if checkModellinTime == True:
	p = subprocess.Popen(['sacct', '-S %s -u jmendieta --format=jobid,state,elapsed' %(date)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	p = subprocess.check_output("sacct -S %s -u jmendieta --format=jobid,state,elapsed" %(date), shell=True)
