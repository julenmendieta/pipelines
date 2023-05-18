#!/bin/bash
# -*- ENCODING: UTF-8 -*-

##===============================================================================
## SLURM VARIABLES
#SBATCH --job-name=Demulti
#SBATCH --cpus-per-task=16
#SBATCH --mem=30G
#SBATCH --time=08:30:00
#SBATCH -p short
#SBATCH -o /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.out  
#SBATCH -e /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.err 

# HOW TO RUN ME
# Update the sampleSheet and path variables. Then run me
# sbatch /home/jmendietaes/programas/pipelines/demulti/bcl2fastq/00_bcl2fastq.sh


module load bcl2fastq2/2.20.0-foss-2018b


sampleSheet=/home/jmendietaes/programas/pipelines/demulti/samples/SampleSheet-run.csv
runfolderDir=/datos/intercambio/DavidL_Amaia/211021_VH00461_74_AAAJFWVHV
outDir=/home/jmendietaes/data/2021/CRISPR/sequencedData/Run74.laura

bcl2fastq -r 10 \
	-p 20 \
	--no-lane-splitting \
	--minimum-trimmed-read-length 35 \
	--mask-short-adapter-read 22 \
	--create-fastq-for-index-reads \
    --sample-sheet ${sampleSheet} \
    --runfolder-dir ${runfolderDir} \
    --output-dir ${outDir} \
    --processing-threads ${SLURM_CPUS_PER_TASK}