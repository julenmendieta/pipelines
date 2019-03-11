import os
from pytadbit import load_structuralmodels
from copy import deepcopy
from pytadbit.modelling.structuralmodels import StructuralModels
from warnings import warn
import sys

def _extend_models(self, models, keep=0):
        """
        add new models to structural models
        """

	new_models = {}
        #nall  = len(self) + len(self._bad_models)
        #self.define_best_models(nall)
        ids = set(self[m]['rand_init'] for m in range(len(self)))
	# create variable for models with same seed
        skip = []
        if len(models) == 0:
            print '##### No models in one case!!!! #####'
        for m in range(len(models)):
	    # we create the list here to retrieve error if there are no models
            #new_models = {}
	    if models[m]['rand_init'] in ids:
                warn('WARNING: found model with same random seed number, '
                     'SKIPPPING')
                skip.append(m)
		print 'skipped'
        inter = [mo for mo in self] + [me for i, me in enumerate(models) if i not in skip ]
        for i, m in enumerate(sorted(inter,
                                     key=lambda x: x['objfun'])):
            new_models[i] = m
        # combine all models
        self = StructuralModels(
            nloci=self.nloci, models=new_models, bad_models=self._bad_models,
            resolution=self.resolution, original_data=self._original_data,
            clusters=self.clusters, config=self._config, zscores=self._zscores,
            zeros=self._zeros, restraints=self._restraints,
            description=self.description)
	# Add new indexes
        for i, m in enumerate(self):
            m['index'] = i
        # keep the same number of best models
        if keep != 0:
            self.define_best_models(keep)
        return self

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
    if len(sys.argv) == 2:  # position zero is extra
        fi = 'lala'
        while fi == 'lala':
            files1 = os.listdir(pathIn)
            for fi1 in files1:
                if fi1.startswith('reg'):
                    fi = fi1
        #fi = '_'.join(fi.split('_')[0].split('Scaled01'))
	flag = '%s_%s' %(fi.split('Scaled01')[0], fi.split('Scaled01')[1].split('_')[0])
        #flag = '%s_%s' %(pathIn.split('/')[-5], fi)
    else:
        flag = sys.argv[2]

print flag
# First store all models in one list
files = os.listdir(pathIn)
modelos = []
for fi in files:
    if fi[-6:] == 'models':
        a = load_structuralmodels(pathIn + fi)
        modelos.append(a)

# then merge them
mods1 = modelos[0]
for i, mods2 in enumerate([m for m in modelos[1:]]):
    #mods1 = _extend_models(mods1, mods2, i + 2)
    mods1 = _extend_models(mods1, mods2)

print len(mods1)
mods1.save_models(pathOut+'%s.modelsAll'%(flag))

