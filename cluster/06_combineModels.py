import os
from pytadbit import load_structuralmodels
from copy import deepcopy
from pytadbit.modelling.structuralmodels import StructuralModels
from warnings import warn
import sys
import subprocess

def _extend_models(models0, models, keep=0):
    """
    add new models to structural models
    """

    new_models = {}
    #nall  = len(models0) + len(models0._bad_models)
    #models0.define_best_models(nall)
    ids = set(models0[m]['rand_init'] for m in range(len(models0)))
    # create variable for models with same seed
    skip = []
    if len(models) == 0:
        print '##### No models in one case!!!! #####'
    # check for duplicated rand_init ids
    for m in range(len(models)):
        # we create the list here to retrieve error if there are no models
        # new_models = {}
        if models[m]['rand_init'] in ids:
            warn('WARNING: found model with same random seed number, '
                    'SKIPPPING')
            ## WITH LAMMPS WE HAVE ANOTHER PARAMETER FOR RANDOM INIT SO NO 
            #SENSE ON REMOVING IN HERE
            skip.append(m)
            print 'Model skipped due to same rand_init'
    inter = [mo for mo in models0] + [me for i, me in 
                                        enumerate(models) if i not in skip]
    for i, m in enumerate(sorted(inter, key=lambda x: x['objfun'])):
        new_models[i] = m
    # combine all models
    models0 = StructuralModels(
        nloci=models0.nloci, models=new_models, bad_models=models0._bad_models,
        resolution=models0.resolution, original_data=models0._original_data,
        clusters=models0.clusters, config=models0._config, zscores=models0._zscores,
        zeros=models0._zeros, restraints=models0._restraints,
        description=models0.description)
    # Add new indexes
    for i, m in enumerate(models0):
        m['index'] = i
    # keep the same number of best models
    if keep != 0:
        models0.define_best_models(keep)
    return models0

##################################################
if len(sys.argv) == 0:
    pathIn = '/scratch/devel/jmendieta/Ferrer/modelling/TBdyn/reg20Norm/scale01/finalModel/'
    pathOut = pathIn
    flag = 'TB_reg20C4.0L-0.5U0.6M300'
#################################################
else:
    pathIn = sys.argv[1]
    pathOut = pathIn
    # build flag
    if len(sys.argv) == 3:
        flag = sys.argv[2]

print flag

# First check the files we have
files = os.listdir(pathIn)
# then check if we have different conformation parameters and split
diffMods = set()
for fi in files:
    if fi[-6:] == 'models':
        diffMods.add('_'.join(fi.split('_')[:-1]))

for di in diffMods:
    # Store all models in one list
    modelos = []
    fileNames = []
    for fi in files:
        if fi[-6:] == 'models' and fi.startswith(di):
            a = load_structuralmodels(pathIn + fi)
            fileNames.append(fi)
            modelos.append(a)

    # get flag
    if len(sys.argv) == 2:  # position zero is extra
        fi = fileNames[0]
        try:
            flag = '%s_%s_%s' %(pathIn.split('/')[-4], 
                        fi.split('_')[0], 
                        fi.split('_')[1].split('Scaled01')[1])
        # If we do have it
        except:
            flag = '%s_%s_%s' %(pathIn.split('/')[-4],
                                            fi.split('_')[0],
                                            fi.split('_')[2].split('Scaled01')[1])
            #flag = '%s_%s' %(pathIn.split('/')[-5], fi)

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
                #mods1 = _extend_models(mods1, mods2, i + 2)
                mods1 = _extend_models(mods1, mods2)
                if len(mods2) == 0:
                    print fileNames[i+ i0 + 1]

    print len(mods1)
    mods1.save_models(pathOut+'%s.modelsAll'%(flag))

# remove all files that are in the merge
clean = True
if clean == True:
    for fi in files:
        if fi[-6:] == 'models':
            p = subprocess.Popen(['rm', '%s%s' %(pathIn, fi)],
                                             stdout=subprocess.PIPE, 
                                             stderr=subprocess.PIPE)
            out, err = p.communicate()
            if len(err) != 0:
                print fi, err
