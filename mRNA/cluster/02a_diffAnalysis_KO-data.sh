#!/bin/bash
# -*- ENCODING: UTF-8 -*-

##===============================================================================
## SLURM VARIABLES
#SBATCH --job-name=mRNA_diff-ko
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --time=00-05:00:00
#SBATCH -p short
#SBATCH -o /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.out  
#SBATCH -e /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.err 
##SBATCH --dependency=afterany:793158



# HOW TO RUN ME
# sbatch /home/jmendietaes/programas/pipelines/mRNA/cluster/02a_diffAnalysis_KO-data.sh \
#/home/jmendietaes/data/2021/mRNA/allProcessed \
#01_projectRestart



##==============================================================================
## GLOBAL VARIABLES

# path where we have the folder structure for chip analysis 
# (inside we have ${basePath}/bamfiles/valid/)
basePath=$1
#basePath=/home/jmendietaes/data/2021/mRNA/allProcessed

## Project name
# We will find salmon output and STAR BAM files in there using our defined
#   folder structure
projectName=$2
#projectName=01_projectRestart

outpath=${basePath}"/furtherAnalysis/${projectName}"
salmonOut=${outpath}/counts
# Path to nextflow scripts
subScripts="/home/jmendietaes/programas/pipelines/mRNA/cluster/sub-scripts"

# Coma separated string with all posible control IDs 
# Used for batch corrected analysis of same KOs
# Set to "no" to not do it
posibleControls="NTC,WT,NTC0005,NtC5,V12h"


##==============================================================================
## Required Software

# Path to R
R="/home/jmendietaes/programas/miniconda3/envs/Renv/bin/Rscript"
# Variables from the job
nCPU=$SLURM_CPUS_PER_TASK

##==============================================================================

# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

##==============================================================================

##############################
# Diff expression with DESeq2
##############################
echo -e "Starting Diff analysis ---------------------------------------------\n"

mkdir -p ${outpath}/DESeq2_batchCorrect 

${R} ${subScripts}/02a_NR_DESeq2_KO-diffExpr.r \
                    --countTable ${salmonOut}/tximportMerge.gene_counts.tsv \
                    --outdir ${outpath}/DESeq2_batchCorrect/ \
                    --outprefix "RNA_DESeq2"\
                    --cores ${nCPU} \
                    --controls ${posibleControls}

mkdir -p ${outpath}/DESeq2_batchCorrect/gatheredDESeq

cd ${outpath}/DESeq2_batchCorrect/gatheredDESeq
ln -s ../*results* . ; 



echo -e "END ------------------------------------------------------------------"

seff $SLURM_JOBID

exit 0
