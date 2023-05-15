#!/bin/bash
# -*- ENCODING: UTF-8 -*-

##===============================================================================
## SLURM VARIABLES
#SBATCH --job-name=tadbitTools
#SBATCH --cpus-per-task=8
#SBATCH --mem=15G
#SBATCH --time=00-48:00:00
#SBATCH -p medium
#SBATCH -o /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.out  
#SBATCH -e /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.err 
#SBATCH --dependency=afterany:802941


##SBATCH --mail-type=END
##SBATCH --mail-user=user@mail.es
# HOW TO RUN ME
# for i in *fastq.gz; do echo $i | sed 's/_R._001.fastq.gz//g' ; done | sort | uniq > samplesNames.txt  
# N=`cat samplesNames.txt | wc -l`
# sbatch --array=1-${N} /home/jmendietaes/programas/pipelines/HiC/cluster/01_tadbitTools.sh

# TADbit tools tutorial (not published)
#https://docs.google.com/document/d/16YseE-v9NwYXZRKKevkGohtJiNdl3VX42wfBPWfkZ-M/edit


# How to run me
# realpath *fastq.gz > toMapp.txt in the location with the fastqs
#and assign it in "filesPath" variable
# python /home/jmendietaes/programas/pipelines/HiC/cluster/00_createMappChia.py
######################  TO CHANGE #####################
# path to file with read1 and read 2 of each experiment
filesPath='/home/jmendietaes/data/2021/HiC/sequencedData/finalDatasets'
#realpath *fastq.gz > toMapp.txt

# load genome path
genome='/home/jmendietaes/referenceGenomes/mm10_reordered/mm10.reordered.fa'
#genome='/home/jmendietaes/referenceGenomes/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa'
#genome='/ssd/genomes/dm6/dm6.fa'

# Variable with chromosome names to analyze
chromCheck=$(for n in {{1..19},X}; do echo chr${n}; done)

# load GEM index path
# GEM2
gem_index_path='/home/jmendietaes/referenceGenomes/mm10_reordered/gem2Index/mm10.reordered.gem'
#gem_index_path='/home/jmendietaes/referenceGenomes/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.gem'
# GEM3
#gem_index_path='/home/jmendietaes/referenceGenomes/mm10_reordered/gem3Index/mm10.reordered.gem'


# Mapper binary
# Gem 2
mapBin="/home/jmendietaes/programas/GEM/GEM-binaries-Linux-x86_64-core_i3-20121106-022124/gem-mapper"
# Gem 3
#mapBin="/home/jmendietaes/programas/gem3-mapper/bin/gem-mapper"
# get output folder
path='/home/jmendietaes/data/2021/HiC/allProcessed/'


# Restriction Enzymes
# None to use iterative mapping, or a string with each of the used enzymes
# separated by "---"
#rEnz = None
rEnz="DpnII DdeI"

##===============================================================================

# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

##===============================================================================


#####################################################
# get number of CPU to use
nthreads=$SLURM_CPUS_PER_TASK

FILES=($(cat ${filesPath}/samplesNames.txt))
filename=${FILES[$SLURM_ARRAY_TASK_ID - 1]}

# get some paths
read1_path="${filesPath}/${filename}_R1_001.fastq.gz"
read2_path="${filesPath}/${filename}_R2_001.fastq.gz"

finalOut=${path}/${filename}
mkdir -p ${finalOut}

echo "Align read 1"
tadbit map -w ${finalOut} --mapper gem --mapper_binary ${mapBin} \
           --fastq ${read1_path} --index ${gem_index_path} --renz ${rEnz} \
           --read 1 -C ${nthreads} --noX  #--mapper_param "--alignment-local-min-identity 15"
echo
echo "Align read 2"
tadbit map -w ${finalOut} --mapper gem --mapper_binary ${mapBin} \
           --fastq ${read2_path} --index ${gem_index_path} --renz ${rEnz} \
           --read 2 -C ${nthreads} --noX  #--mapper_param "--alignment-local-min-identity 15"
echo
echo "Parsing"
tadbit parse -w ${finalOut} --genome ${genome} --noX --compress_input
echo
echo "Filtering"
tadbit filter -w ${finalOut} -C ${nthreads} --noX  --apply 1 2 3 4 7 9 10 \
            --clean --valid --compress_input
echo
echo "normalize"
tadbit normalize -w ${finalOut} --resolution 50000 -C ${nthreads} --noX --valid
tadbit normalize -w ${finalOut} --resolution 100000 -C ${nthreads} --noX --valid
tadbit normalize -w ${finalOut} --resolution 10000 -C ${nthreads} --noX --valid

# Call TADs at 50kb. Use TADbit alg to call TADs, combine all the borders in all
# conditions, and then look at the Insulation Score (SC) there
# calculas el insulation score sobre todo el genoma,luego haces un par 
#de piruetas y predices borders... pero los borders de tadbit son mejores. 
#por eso lo que yo hago es usar los borders de tadbit, y luego si necesito 
#ver valores de insulacion continuos, uso el insulation score
# Get insulation score in pytadbit
#http://3dgenomes.github.io/TADbit/tutorial/tutorial_8-Compartments_and_TADs_detection.html#insulation-score
# echo
echo "TAD calling"
for chrom in ${chromCheck}; do 
    tadbit segment -w ${finalOut} -r 50000 --noX --only_tads -C ${nthreads} -c $chrom; 
done
# Call Compartments at 100kb 
echo
echo "Compartments"
for chrom in ${chromCheck}; do
    tadbit segment -w ${finalOut} -r 100000 --noX --only_compartments \
                -C ${nthreads} --ev_index 1 --smoothing_window 3 \
                -c $chrom --savecorr --fasta ${genome}
    tadbit segment -w ${finalOut} -r 100000 --noX --only_compartments \
                -C ${nthreads} --ev_index 2 --smoothing_window 3 \
                -c $chrom --savecorr --fasta ${genome}
done

# From here we go to pyhton
#http://3dgenomes.github.io/TADbit/tutorial/tutorial_8-Compartments_and_TADs_detection.html#insulation-score

# Take a look
#tadbit describe -w PATH -t job
#tadbit bin -w ${finalOut} -r 10_000 -c chr1:79_600_000-80_600_000 --interactive --triangular --only_plot --tad_def 54 --norm 'norm'
# Copy output 
cat /home/jmendietaes/jobsSlurm/outErr/tadbitTools_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out >> ${finalOut}/${filename}_TADbit_output.txt

