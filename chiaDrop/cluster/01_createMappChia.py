import time

######################  TO CHANGE #####################
# path to file with read1 and read 2 of each experiment
filesPath = '/home/jmendietaes/data/2021/chia-drop/fastq/toMapp.txt'

# load genome path
genome = '/home/jmendietaes/referenceGenomes/mm10_reordered/mm10.reordered.fa'
#genome = '/home/jmendietaes/referenceGenomes/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa'
#genome = '/ssd/genomes/dm6/dm6.fa'

# load GEM index path
gem_index_path = '/home/jmendietaes/referenceGenomes/mm10_reordered/mm10.reordered.gem'
#gem_index_path= '/home/jmendietaes/referenceGenomes/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.gem'

# get output folder
path = '/home/jmendietaes/data/2021/chia-drop/mapOut/'

# get temporal output folder if present
tempOut = '/datos/mapTemp/'
#tempOut = '/home/jmendietaes/data/2021/chia-drop/mapOut/'

# minimum allowed number of bp in a range
minseq = 25

# maximum allowed number of bp to be mapped (max length of read=)
maxseq = 130

# get number of CPU to use
nthreads = 16

# Path for location of scripts
scriptsPath = ''

#####################################################

reads = []
with open(filesPath, 'r') as f:
    for line in f:
        read1 = line.split()
        read2 = f.next().split()
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
parser.add_argument('-fq1','--fastq1', help = 'Path to fasq file of read 1',required=True)
parser.add_argument('-fq2','--fastq2', help = 'Path to fasq file of read 2', required=True)
parser.add_argument('-gnm','--genome', help = 'Path to reference genome', required=True)
parser.add_argument('-gi','--gemindex', help = 'Path to GEM index file', required=True)
parser.add_argument('-op','--outpath',help='Path to output folder', required=True)
parser.add_argument('-opt','--tempdir',help='Path for temporal output folder', required=False)
parser.add_argument('-mins','--minseq',help='Minimum bp length to map (def. 25bp)', required=False)
parser.add_argument('-maxs','--maxseq',help='Maximum bp length to map (Min read length)', required=True)
parser.add_argument('-t','--threads',help='Number of CPU to use', required=False)

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
#SBATCH --cpus-per-task=16
#SBATCH --mem=10G
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

python $orden''' %(njobs, runfile)

fout.write(cmd)
