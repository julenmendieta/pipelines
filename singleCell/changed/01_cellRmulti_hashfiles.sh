#!/bin/bash
# -*- ENCODING: UTF-8 -*-

##===============================================================================
## SLURM VARIABLES
#SBATCH --job-name=cellRmulti
#SBATCH --cpus-per-task=16
#SBATCH --mem=40Gb
#SBATCH --time=20:00:00
#SBATCH -p short
#SBATCH -o /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.out  
#SBATCH -e /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.err 
##SBATCH --dependency=afterany:599711


# original ram was 258Gb, and cpu 24

# HOW TO RUN ME
#sbatch /home/jmendietaes/programas/PhD/singleCell/changed/01_cellRmulti_hashfiles.sh

# OBJECTIVE
# run cellRanger to analyze cell hashing data
# https://kb.10xgenomics.com/hc/en-us/articles/360028612751-Can-I-use-Cell-Ranger-to-analyze-cell-hashing-data-

#source $CODEBASE/tfcf/setup.sh

#### TO CHANGE
# Main data project file
basepath="/home/jmendietaes/data/2021/singleCell"


# base name of the files to check (separated by space)
# avoid adding _Library.csv, it will be added later
filesCheck="SAOS_unperturbed"

# Genome ID
#genomeId="refdata-gex-mm10-2020-A"
genomeId="refdata-gex-GRCh38-2020-A"

# Path where we will store output data
outputPath="${basepath}/allProcessed"

# Path where we have the resource files from cellRanger
genomeBase=$basepath/additionalFiles/omicstmp/
rangerPath="${basepath}/allProcessed/rangerFiles"
# path for scripts location
scriptsPath="/home/jmendietaes/programas/PhD"

featurePath=/home/jmendietaes/data/2021/singleCell/allProcessed/rangerFiles/SAOS_unperturbed_reference.csv

#### MAIN CODE
# load modules
module load CellRanger/6.1.1

echo ${SLURM_MEM_PER_NODE}
adjustedMem=$(echo "print(int(round($SLURM_MEM_PER_NODE*0.98,0)/1000))" | python3)
# fast check in cases where RAM given in Gb (value will be lower than 1 in most cases)
# 100Gb / 1000 = 0.1
if [ "$adjustedMem" -lt "5" ]; then
    adjustedMem=$(echo "print(int(round($SLURM_MEM_PER_NODE*0.98,0)))" | python3)
fi  

# If genome doesn't exist, then create it
if [ ! -f "${genomePath}/genome.fa" ]; then
    bash ${scriptsPath}/singleCell/changed/00_PREP_Genome.sh
fi

if [ ! -f "${outputPath}/omicstmp" ]; then
    mkdir -p ${outputPath}/omicstmp
fi


########## BASIC ANALYSIS
genomePath="${genomeBase}/${genomeId}"
cd $outputPath/omicstmp

for id in ${filesCheck}; do

    echo $id
    
    # check that feature file exist
    if [ ! -e ${featurePath} ]; then
        echo "You forgot to create the features files"
        exit 1
    fi


    echo "Dataset"
    #cat $rangerPath/${id}_dataset.csv

    # cellranger multi --id=$id --csv=$rangerPath/${id}_dataset.csv \
    #     --localcores=${SLURM_CPUS_PER_TASK} \
    #     --localmem=${adjustedMem} &> ${id}.log


    cellranger count --id=$id \
     --no-bam \
     --libraries=$rangerPath/${id}_Library.csv \
     --transcriptome="${genomePath}-Extended/" \
     --feature-ref=${featurePath} \
     --localcores=${SLURM_CPUS_PER_TASK} \
     --localmem=${adjustedMem} \
     --expect-cells=10000 &> ${id}.log


    mkdir -p ${outputPath}/Data/$id
    mv $id/outs ${outputPath}/Data/$id
done


