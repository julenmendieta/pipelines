#!/bin/bash
# -*- ENCODING: UTF-8 -*-

##===============================================================================
## SLURM VARIABLES
#SBATCH --job-name=cellRaggr
#SBATCH --cpus-per-task=24
#SBATCH --mem=30Gb
#SBATCH --time=10:00:00
#SBATCH -p short
#SBATCH -o /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.out  
#SBATCH -e /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.err 
##SBATCH --dependency=afterany:599711


# original ram was 258Gb, and cpu 24

# HOW TO RUN ME
#sbatch /home/jmendietaes/programas/PhD/singleCell/changed/02_cellRaggregate.sh

# OBJECTIVE
# aggregate cellranger count experiments

#### TO CHANGE
# Folder were we have the aggregate csv files
#https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/aggregate
inpath="/home/jmendietaes/data/2021/singleCell/allProcessed/aggregate/csv_files"

# Path where we will store output data
outputPath="/home/jmendietaes/data/2021/singleCell/allProcessed/aggregate/aggregates"

# name of the files to check (separated by space). ID will be the file name minus .csv
filesCheck="inVivo_14d.csv inVivo_28d.csv"


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


########## BASIC ANALYSIS
cd $outputPath

for file1 in ${filesCheck}; do
    id=${file1::-4}
    echo $id

    cellranger aggr --id=$id \
     --csv=${inpath}/${file1} \
     --normalize=mapped \
     --localcores=${SLURM_CPUS_PER_TASK} \
     --localmem=${adjustedMem} \
      &> ${outputPath}/${id}.log

    
done


