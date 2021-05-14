import os
from copy import deepcopy
from warnings import warn
import sys
import subprocess
from taddyn.utils.modelAnalysis     import save_models
import copy
import cPickle

def extend_models(modelsBase, modelsAdd, keep=0):
    """
    combine taddyn models
    """
    
    ## Copy base models dictionary
    modelsFinal = copy.deepcopy(modelsBase)

    ## get all rand init from previous models
    randInit_indexes = []
    for mo in modelsBase['models']:
        randInit_indexes += [modelsBase['models'][mo]['rand_init']]

    ## now add models if dont share rand_init
    index = len(randInit_indexes)
    duplicated = 0
    for mo in modelsAdd['models']:
        randinitAdd = modelsAdd['models'][mo]['rand_init']
        if randinitAdd not in randInit_indexes:
            modelsFinal['models'][index] = modelsAdd['models'][mo]
            index += 1
        else:
            print('Model skipped due to same rand_init: %s' %(randinitAdd))
            duplicated += 1
    
    if duplicated != 0:
        print('%s models were skipped due to repteated rand_init' %(duplicated))
    
    return modelsFinal


## Get input path
#parser = argparse.ArgumentParser(description='')
#parser.add_argument('-p','--path', help='modelsPath',required=True)
#args = parser.parse_args()
#pathIn=float(args.path)
pathIn = sys.argv[1]
pathOut = pathIn

## Get different models in folder and merge the ones computed
## with same parameters
# First check the files we have
files = os.listdir(pathIn)
# then check if we have different conformation parameters and split
diffMods = set()
for fi in files:
    if fi.endswith('modelsTemp'):
        diffMods.add('_'.join(fi.split('_')[:-1]))


for di in diffMods:
    # Store all models in one list
    modelos = []
    fileNames = []
    for fi in files:
        # here two things can happen;
        # i) all models are from direct output
        # ii) i changed the name of the merge output from
        #...modelsAll to ...models to merge new build models
        # case i)
        if fi.endswith('modelsTemp') and fi.startswith(di):
            try:
                with open(pathIn + fi, "rb") as input_file:
                    a = cPickle.load(input_file)
                fileNames.append(fi)
                filenameForFlag = fi
                modelos.append(a)
            except:
                print('Empty file:')
                print(pathIn + fi)
            
        # nnnnnn maybe this change to modelsAll
        # If we have previous finished models we want to combine them
        if fi.endswith('models') and fi.startswith(di):
            print('previous model ensemble found, combining')
            print(pathIn + fi)
            with open(pathIn + fi, "rb") as input_file:
                a = cPickle.load(input_file)
            fileNames.append(fi)
            modelos.append(a)
    
    # get flag
    fi = filenameForFlag
    #rint fi
    flag = '%s_%s_%s' %(fi.split('_')[0], 
        fi.split('_')[1], 
        fi.split('_')[2])

    print(flag)
    ## then merge them
    # check that first model has data
    full = False
    stop = False
    i0 = 0
    while full == False and stop == False:
        mods1 = modelos[i0]
        if len(mods1) == 0:
            print '##### No models in one case!!!! #####'
            print fileNames[i0]
            i0 += 1
        else:
            full = True
        # if reach end of models and all have no data
        if i0 >= len(modelos):
            stop = True
    # if we have data
    if full == True:
        # can be that just have one model
        if i0 == len(modelos) - 1:
            mods1 = modelos[i0]
        # if we have more than one
        else:
            for i, mods2 in enumerate([m for m in modelos[i0 + 1:]]):
                mods1 = extend_models(mods1, mods2)
                if len(mods2) == 0:
                    print('Empty model reached last combination phase')
                    print(fileNames[i+ i0 + 1])

    # remove combined files before storing new
    clean = True
    if clean == True:
        print('--- Deleting .model files ---')
        for fi in fileNames:
            print 'Delete: %s' %fi
            p = subprocess.Popen(['rm', '%s%s' %(pathIn, fi)],
                                           stdout=subprocess.PIPE, 
                                           stderr=subprocess.PIPE)
            out, err = p.communicate()
            if len(err) != 0:
               print fi, err

    print('%s models joined' %len(mods1['models']))
    save_models(mods1, pathOut+'%s.models'%(flag))
