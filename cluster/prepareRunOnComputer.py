import multiprocessing 
import os
from time import sleep

loadMatrixFile = '/scratch/julen/PCHiC_assesment/modelling/matrixList.txt'
n_cpus = 8
totalProcesors = 32

maxJobs = (totalProcesors / n_cpus) - 1

# Get matrices from file
matPaths= []
with open(loadMatrixFile,'r') as f:
    for line in f:
        matPaths += line.split()

commands = []
for matPath in matPaths:
    path = '/'.join(matPath.split('/')[:-1]) + '/'
    # get file with running istructions
    files = os.listdir(path)
    for fi in files:
        if fi.endswith('array'):
            runFile = fi
    runFile = path + runFile

    # get commands
    with open(runFile, 'r') as f:
        for line in f:
            commands.append(line.rstrip())

# Create queue
the_queue = multiprocessing.Queue()

# Create the function that will run the jobs
def worker2(queue):
    print os.getpid(),"working"
    while True:
        item = queue.get(True)
        print os.getpid(), "got", item
        os.system(item)
        return item

# Create pool object
the_pool = multiprocessing.Pool(4, worker2,(the_queue,))
#                        don't forget the coma here  ^


# run jobs
for co in commands:
    #print 'python %s' %co
    item = the_queue.put('python %s' %co)

   
# Wail till all jobs have finished
waiting = the_queue.qsize()
while waiting > 0:
    sleep(5)
