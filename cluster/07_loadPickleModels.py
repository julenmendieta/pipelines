from pytadbit.modelling.structuralmodels import StructuralModels
from cPickle import load

# if something fails in further analysis quite likely is because the areas we put as None in StructuralModels
#should have another value

## input
paths = [] # paths for the model pickle files
nloci = 0  # number of loci in our models
resol = 5000  # resolution
originalMatrix = [[]]  # interaction matrix to which the models will be compared. mods.correlate_with_real_data(
pathOut = '.models'  # path and name file were to write models

## run
results = []
# load pickles
for model_path in paths:
    with open(model_path, "rb") as input_file:
        m = load(input_file)
    results.append((m[0], m[1]))

# load all models and sort by objective function
models = {}
for i, (_, m) in enumerate(
    sorted(results, key=lambda x: x[1][0]['objfun'])[:]):
    models[i] = m[0]

# add index and descritption
for i, m in enumerate(models.values()):
    m['index'] = i
    m['description'] = 'cute model'


# create structural models object
models0 = StructuralModels(
        nloci=nloci, models=models, bad_models=None,
        resolution=resol, original_data=originalMatrix,
        clusters=None, config=models0._config, zscores=None,
        zeros=None, restraints=None,
        description=None)

# save models
models0.save_models(pathOut)