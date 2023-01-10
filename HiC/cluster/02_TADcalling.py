from pytadbit.tadbit import insulation_score, insulation_to_borders
from pytadbit.parsers.hic_parser import load_hic_data_from_bam
from pytadbit import Chromosome
import sys
import os
import pickle

# python /home/jmendietaes/programas/PhD/HiC/cluster/02_TADcalling.py

# List of chromosome to analyse
chroms = [f'chr{i}' for i in range(1, 20)] + ['chrX', 'chrY']

# BAM file path
bampath = "/home/jmendietaes/data/2021/HiC/mapOut/ToolsOut/DM_HiC_normal_S34/03_filtered_reads/intersection_94e4921e80.bam"

# Resolution at which to call TADs
resol = 50000

# Normalisation biases path
biasesP = "/home/jmendietaes/data/2021/HiC/mapOut/ToolsOut/DM_HiC_normal_S34/04_normalization/biases_50kb_17d8f87d9d.pickle"

# Path to output directoy
outpath = "/home/jmendietaes/data/2021/HiC/mapOut/ToolsOut/DM_HiC_normal_S34/06_segmentation/tadsPy_50kb"

# ncpus to use at the time to load the matrix and call TADs
ncpus = 4

##########################  Start program  ##################

# Create output directory if it doesnt exist
if not os.path.exists(outpath):
    os.makedirs(outpath)


toStore = {}
for chrom in chroms:
    toStore[chrom] = {}
    
    print(chrom)
    # Open bam file at chrom location
    hic_data = load_hic_data_from_bam(bampath, resol, 
                                        region=chrom, 
                                        filter_exclude=(),
                                        ncpus=ncpus)
    # Normalise by decay
    hic_data.load_biases(biasesP)
    hic_data.normalize_expected()

    # Get insulation score
    # the square size should be 500 kb as close as possible from the diagonal
    wsize = (1, 500000//resol)
    # the delta is to look for increases in insulation around a given bin. 
    # Should be around 100 kb, but more than one bin
    delta = max(100000//resol, 2)
    # Compute insulation score
    insc1, delta1 = insulation_score(hic_data, [wsize], 
                    resolution=resol, normalize=True, delta=delta)

    # Search for borders
    borders1 = insulation_to_borders(insc1[wsize], delta1[wsize], 
                                    min_strength=0.1)

    # Store in variable
    toStore[chrom]['ISscores'] = insc1
    toStore[chrom]['ISdelta'] = delta1
    toStore[chrom]['ISborders'] = borders1

    # Call TADs with TADbit algorithm (already done with TADbit tools)
    # crm = Chromosome(chrom)
    # crm.add_experiment(chrom,
    #                 hic_data=[hic_data.get_matrix(focus=chrom)],
    #                 norm_data=[hic_data.get_matrix(focus=chrom,
    #                                             normalized=True)],
    #                 resolution=resol)
    # crm.find_tad([chrom], n_cpus=ncpus)

    # # Now we store TAD borders
    # crm.experiments[0].write_tad_borders(density=True, 
    #                         savedata=f"{outpath}/{chrom}-{resol}_TADbit-TADborders.tsv", 
    #                         normalized=True)


with open(f"{outpath}/InsulationOut_{resol}.pickle", 'wb') as handle:
    pickle.dump(toStore, handle, protocol=pickle.HIGHEST_PROTOCOL)

#with open(f"{outpath}/InsulationOut_{resol}.pickle", 'rb') as handle:
#    toStore = pickle.load(handle)