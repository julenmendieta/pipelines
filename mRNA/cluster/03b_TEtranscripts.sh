#!/bin/bash
# -*- ENCODING: UTF-8 -*-

##==============================================================================
## SLURM VARIABLES
#SBATCH --job-name=TEtranscripts
#SBATCH --cpus-per-task=1
#SBATCH --mem=30G
#SBATCH --time=00-24:00:00
#SBATCH -p short
#SBATCH -o /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.out  
#SBATCH -e /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.err 
##SBATCH --dependency=afterany:793158


##SBATCH --mail-type=END
##SBATCH --mail-user=user@mail.es
# HOW TO RUN ME
# sbatch /home/jmendietaes/programas/pipelines/mRNA/cluster/03b_TEtranscripts.sh \
#/home/jmendietaes/data/2021/mRNA/allProcessed/TEtranscripts \
#/home/jmendietaes/referenceGenomes/mm10_reordered/STAR/mm10.reordered \
#03_TEtranscripts
#NTC

# OBJECTIVE
# Count transcripts on TE families

# NOTES
# "In our experience, we recommend around 20-30Gb of memory for analyzing 
#   human samples (hg19) with around 20-30 million mapped reads when running 
#   on a cluster."
# Code based on notes from:
# https://doi.org/10.1093/bioinformatics/btv422
# https://doi.org/10.1007/978-1-4939-7710-9_11

##==============================================================================
## GLOBAL VARIABLES


## path to the main fully processed data output
basePath=$1
#basePath="/home/jmendietaes/data/2021/mRNA/allProcessed/TEtranscripts"

##  REFERENCE Genome info
REFERENCE_DIR=$2
#REFERENCE_DIR="/home/jmendietaes/referenceGenomes/mm10_reordered/STAR/mm10.reordered"
# REFERENCE_DIR includes the path and naming base of the reference genome, we will
#   add .sizes, .blacklist.bed to the base to get the rest of paths. Genome
#   indexes must follow the same base name

GenomeIndex=$REFERENCE_DIR

# Get all kind of posible control IDs (comma separated)
posibleControls=$4
#posibleControls="NTC,WT,NTC0005,NTC5,V12h"


## Extra path varaibles
# GTF to TE
TEgtf=/home/jmendietaes/referenceGenomes/annotRegions/mm10/mm10_rmsk_TE.gtf

# Define GTF to make STAR index (at $REFERENCE_DIR path in genes folder)
indexOutP=$(dirname ${GenomeIndex})
useGTF=$(realpath ${indexOutP}/genes/*gtf)
# use singleCell GTF for test
#useGTF=/home/jmendietaes/data/2021/singleCell/additionalFiles/refdata-gex-mm10-2020-A/genes/genes_mainGenes.gtf

## Project name
# We will find salmon output and STAR BAM files in there using our defined
#   folder structure
projectName=$3
#projectName=03_TEtranscripts

# Were to look for bam files
bamsPath="${basePath}/bamfiles/wholeGenome/valid/${projectName}"
outpath=${basePath}"/furtherAnalysis/${projectName}"


##==============================================================================
## Required Software
module load SAMtools/1.12-GCC-10.2.0
# GCC needed for STAR
module load GCCcore/11.2.0
# this is for bamCoverage, but i already have it in conda and newer version
#module load deepTools/3.2.0-foss-2018b-Python-2.7.15


# Path to STAR binary
TEtranscriptsP="/home/jmendietaes/programas/miniconda3/bin/TEtranscripts"
TEcountP="/home/jmendietaes/programas/miniconda3/bin/TEcount"
# Export R path
export PATH="/home/jmendietaes/programas/miniconda3/envs/Renv/bin:$PATH"


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

# function to check if the given first file doesnt exist or is older than 
# the second input file
fileNotExistOrOlder () {
    # check if the file exists of it was created with a previous bam version 
    analyse="no"
    if [ ! -e $1 ]; then
        analyse="yes"
    # only proceed if the output file is older than the bam file
    # in this way if we resequenced and kept the name the analysis 
    # will be repeated
    else
        for tfile in $2; do
            # If $1 is older than any $2
            if [[ $1 -ot ${tfile} ]] ; then
                analyse="yes"
                echo $1" older than "${tfile}
            fi
        done
    fi
}

##==============================================================================

if [ ! -e ${outpath}/TEtranscript ]; then
    mkdir -p ${outpath}/TEtranscript
fi



#################################
## Define strandness of samples
#################################

samples=$(for i in $bamsPath/*bam; do echo $i | sed "s+${bamsPath}\/++g" | \
                                       sed 's/.STAR.sort.rmchr.bam//g'; done)

lTypes=$(for s in ${samples}; do 
    grep -A 1 "LIBRARY TYPE" ${basePath}/QC/allStepsSummary/summary_${s}_TE.txt |\
    tail -n 1; done | sort | uniq)
nTypes=$(echo ${lTypes} | wc -w)
if [ "$nTypes" -gt 1 ] ; then
    echo "ERROR: More than one library type"
    echo $nTypes 
    exit 1
fi
# Define strandness
# 0 (unstranded), 1 (stranded) and 2 (reversely stranded)
lType_="no"
if [[ ${lTypes} == *"SF" ]] ; then
    lType_="forward"
elif [[ ${lTypes} == *"SR" ]] ; then
    lType_="reverse"
fi




#################################
## TEtranscriptsP: TE families that are enriched
#################################

# Replicate samples of the same condition should be input together as a 
#   group and separated by a space:
# TEtranscripts --format BAM --stranded reverse -t piwiKD_1_ Aligned.out.bam \
#  piwiKD_2_Aligned.out.bam -c control_1_Aligned.out.bam \
#  control_2_Aligned.out.bam --GTF ~/Index/gene_ann.gtf \
#  --TE ~/Index/TE_ann.gtf --mode multi --project piwiKD_vs_control \
#  --minread 1 -i 10 --padj 0.05
# TEtranscripts returns two output tables from DESeq: the standard output 
#  table of fold change and p-value statistics for all genes and TEs, 
#  as well as a second table of only those genes and TEs calculated to 
#  be statistically significant in their differential expression between 
#  conditions
# Output:
# piwiKD_vs_control.cntTable contains the estimated raw abundance counts for 
#  all genes and TEs as a tab-delimited table. Each row represents a gene or TE, 
#  each column is a sample, and each value is a raw count.
# piwiKD_vs_control_gene_TE_analysis.txt contains the differential expression 
#  results from DESeq for all genes and TEs.
# piwiKD_vs_control_sigdiff_gene_TE.txt contains a subset of the differential 
#  expression analysis table for only those genes and TEs that passed the 
#  P-value significance criteria

# Comparison is done by cells
allbams=$(find ${bamsPath}/*bam -printf "${bamsPath}/%f\n")
allCells=$(\
    for filename in ${allbams}; do 
        filename=$(basename ${filename})
        mapLib=(${filename//_/ }); 
        mapLib=${mapLib[0]}; 
        mapLib=(${mapLib//-/ }); 
        echo ${mapLib[0]}; done | sort | uniq)

for cell in ${allCells}; do
   cellBams=$(echo $allbams | tr ' ' '\n' | \
            { grep -e "${cell}-\|${cell}_" || :; })

   # Get guides
   allGuides=$(\
         for filename in ${cellBams}; do 
            filename=$(basename ${filename})
            mapLib=(${filename//_/ }); 
            mapLib=${mapLib[0]}; 
            mapLib=(${mapLib//-/ }); 
            echo ${mapLib[1]}; done | sort | uniq)

   grepControls=$(echo ${posibleControls} | sed 's/,/\\\|/g')
   controlBams=$(echo "${cellBams}" | grep -e $grepControls)

   koGuides=$(echo $allGuides | tr ' ' '\n' | grep -v -e $grepControls)

   for ko in ${koGuides}; do
      koBams=$(echo "${cellBams}" | grep -e $ko)

      # Check if we have already a log for this check
      fileNotExistOrOlder "${outpath}/TEtranscript/${cell}_${ko}-vs-Control.cntTable" \
                        "${koBams} ${controlBams}"
      if [[ ${analyse} == "yes" ]]; then
         echo "Runing analysis of:"
         echo "KO: ${koBams}"
         echo "Control: ${controlBams}"
         ${TEtranscriptsP} --sortByPos --format BAM --mode multi \
                  -t ${koBams} -c ${controlBams} \
                  --project "${cell}_${ko}-vs-Control" \
                  --GTF ${useGTF} --TE ${TEgtf} \
                  --stranded ${lType_} --outdir ${outpath}/TEtranscript \
                  --padj 0.05 --foldchange 1.5
      else
         echo "Files already analysed"
         echo "KO: ${koBams}"
         echo "Control: ${controlBams}"
      fi

   done
done


# ${TEcountP} --sortByPos --format BAM --mode multi -b ${bam1} \
#             --GTF ${useGTF} --TE ${TEgtf} \
#             --stranded ${lType_} --outdir ${outDir} --project test1Count





echo -e "END ------------------------------------------------------------------"

seff $SLURM_JOBID

exit 0