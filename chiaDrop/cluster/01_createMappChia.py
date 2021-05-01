import time

######################  TO CHANGE #####################
# path to file with read1 and read 2 of each experiment
filesPath = '/home/jmendietaes/data/2021/chia-drop/fastq/toMapp.txt'
#realpath *fastq.gz > toMapp.txt

# load genome path
#genome = '/home/jmendietaes/referenceGenomes/mm10_reordered/mm10.reordered.fa'
genome = '/home/jmendietaes/referenceGenomes/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa'
#genome = '/ssd/genomes/dm6/dm6.fa'

# load GEM index path
#gem_index_path = '/home/jmendietaes/referenceGenomes/mm10_reordered/mm10.reordered.gem'
gem_index_path= '/home/jmendietaes/referenceGenomes/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.gem'

# get output folder
path = '/home/jmendietaes/data/2021/chia-drop/mapOut/'

# get temporal output folder if present
tempOut = '/datos/intercambio/mapTemp/'
#tempOut = '/home/jmendietaes/data/2021/chia-drop/mapOut/'

# minimum allowed number of bp in a range
minseq = 25

# maximum allowed number of bp to be mapped (max length of read=)
maxseq = 130

# get number of CPU to use
nthreads = 16

# Path for location of scripts
scriptsPath = '/home/jmendietaes/programas/PhD/chiaDrop/cluster/'

#####################################################

reads = []
with open(filesPath, 'r') as f:
    for line in f:
        read1 = line.split()[0]
        read2 = f.readline().split()[0]
        reads += [(read1, read2)]



#flag = matPath.split('/')[-2]
#path='/'.join(matPath.split('/')[:-1]) + '/'
#lammpsOut="/scratch_tmp/%s/temp" %(flag)



# Keep record of the number of jobs to run
njobs = len(reads)
print(njobs)

fecha = time.strftime("%d-%m-%Y")
# File were we are going to keep all the parameter combinations
runfile = '%s%s.Mappsarray'%(path, fecha)
## Create the file with the commands to be run in the array
fout=open(runfile,'w')
for nr in range(0, njobs):
	cmd=''
	cmd+='%s01_mapper.py -fq1 %s '%(scriptsPath, reads[nr][0])
	cmd+= '-fq2 %s '%reads[nr][1]
	cmd+= '-gnm %s '%genome
	cmd+= '-gi %s '%gem_index_path
	cmd+= '-op %s '%path
	cmd+= '-opt %s '%tempOut
	cmd+= '-mins %s '%minseq
	cmd+= '-maxs %s '%maxseq
	cmd+= '-t %s'%nthreads
	#cmd+= '-t %s '%jobTime  # will add it in the arrayMod.jobs so we can update it later
	cmd+= '\n'
	fout.write(cmd)
fout.close()



## Create the file to launch the array (sbatch file)
fout=open('%s/arrayMappChia.cmd'%(path),'w')
cmd=''
cmd+='''#!/bin/bash

#SBATCH --time=96:00:00
#SBATCH -o /home/jmendietaes/jobsSlurm/outErr/%%x_%%A_%%a.out
#SBATCH -e /home/jmendietaes/jobsSlurm/outErr/%%x_%%A_%%a.err
#SBATCH --job-name=mappingChIA-drop
#SBATCH -p medium
#SBATCH --nodes=1
#SBATCH --cpus-per-task=%s
#SBATCH --mem=20G
#SBATCH --array=1-%s%%5

module load GLib/2.54.3-GCCcore-7.3.0
module load GCC

# File were we have located our array commands
file=%s

# Get each command from the file and run them with python
orden=`sed "${SLURM_ARRAY_TASK_ID}q;d" $file`
# will add the command for the temporal folder
orden=`echo $orden`
echo ${SLURM_ARRAY_TASK_ID}
echo ${orden}

python $orden''' %(nthreads, njobs, runfile)

fout.write(cmd)
