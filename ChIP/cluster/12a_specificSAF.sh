#!/bin/bash
# -*- ENCODING: UTF-8 -*-

##===============================================================================
## SLURM VARIABLES
#SBATCH --job-name=specificSAF
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH -p short
#SBATCH -o /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.out  
#SBATCH -e /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.err 


# HOW TO RUN ME
#sbatch /home/jmendietaes/programas/pipelines/ChIP/cluster/12a_specificSAF.sh

# OBJECTIVE
# Get all peak files from folder and generate consensus coordinate SAF file

# Main peak folder
# All files ending by Peak will be used for the consensus (we expect MACS like format)
peakPath=/home/jmendietaes/data/2021/chip/allProcessed/furtherAnalysis/08_projectRestart/peakCalling/peakSelection
# path for the location of the pipeline scripts
scriptsPath="/home/jmendietaes/programas/PhD"

##===============================================================================

# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command failed with exit code $?."' EXIT

##===============================================================================

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

##===============================================================================


echo -e "Starting consensus peak analysis -------------------------------------\n"


chip="allmerged"
peaktype=broadPeak

# select the columns for broadPeak (since they are also included in narrow ones)
mergecols=`seq 2 9 | tr '\n' ','`
expandparam=''

# get all non empty peak files (avoid controls)
peakFiles=$(find -L ${peakPath}/peaks/*Peak -maxdepth 1  -type f ! -size 0 | \
            { grep -v -e "_input" -v -e "_IgG" || :; } )
fileLabels=$(for f in $peakFiles; do echo ${f##*/} |sed "s/_peaks.${peaktype}//g"; done | tr '\n' ',')

# check if the file exists or it was created with a previous peaks version 
prefix="${chip}_${peaktype}_consensusPeaks"
fileNotExistOrOlder "${peakPath}/consensus/${prefix}.saf" "${peakFiles}"
# this outputs analyse as yes or no in lowercase
if [[ ${analyse} == "yes" ]]; then

    sort -T '.' -k1,1 -k2,2n ${peakFiles} \
        | mergeBed -c $mergecols -o collapse > ${peakPath}/consensus/${prefix}.txt

    python ${scriptsPath}/ChIP/cluster/02_NR_macs2_merged_expand.py ${peakPath}/consensus/${prefix}.txt \
        ${fileLabels} \
        ${peakPath}/consensus/${prefix}.boolean.txt \
        $expandparam

    consensusPeakBed=${peakPath}/consensus/${prefix}.bed
    awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print $1, $2, $3, $4, "0", "+" }' \
        ${peakPath}/consensus/${prefix}.boolean.txt > ${consensusPeakBed}
    echo -e "GeneID\tChr\tStart\tEnd\tStrand" > ${peakPath}/consensus/${prefix}.saf
    awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print $4, $1, $2, $3,  "+" }' \
        ${peakPath}/consensus/${prefix}.boolean.txt >> ${peakPath}/consensus/${prefix}.saf
    
fi
    
