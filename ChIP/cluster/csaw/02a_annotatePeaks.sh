#!/bin/bash
# -*- ENCODING: UTF-8 -*-

##===============================================================================
## SLURM VARIABLES
#SBATCH --job-name=AnnotatePeaks
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=0-06:00:00
#SBATCH -p short
#SBATCH -o /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.out  
#SBATCH -e /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.err 


# HOW TO RUN ME
#sbatch /home/jmendietaes/programas/PhD/ChIP/cluster/csaw/02a_annotatePeaks.sh 

# OBJECTIVE
# Annotate peakfiles in a folder

basePath="/home/jmendietaes/data/2021/chip/allProcessed/furtherAnalysis/subsampled_noIgG/binCalling"
inpath="${basePath}/binnedPeaks/consensus"
outpath="${basePath}/Annot/consensus"

# path for the location of the pipeline scripts
scriptsPath="/home/jmendietaes/programas/PhD"
# species shortcut for MACS
species="mm"
speciesGenome="mm10"

## load modules
export PATH="/home/jmendietaes/programas/miniconda3/envs/Renv/bin:$PATH"
export PATH="/home/jmendietaes/programas/miniconda3/bin:$PATH"

module load Homer/4.10-foss-2018b

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
            if [[ $1 -ot ${tfile} ]] ; then
                analyse="yes"
                echo $1" older than "${tfile}
            fi
        done
    fi
}

############################################################
# Annotate consensus peaks with HOMER, and add annotation to boolean output file
############################################################

if [ ! -e ${outpath} ]; then
	mkdir -p ${outpath}
fi

echo -e "Starting consensus peak annotations ------------------------------\n"

consensusFiles=$(find ${inpath}/*tsv -printf "${inpath}/%f\n" | \
            tr '\n' ' ')
chip="allmerged"
for binnedPeaks in ${consensusFiles}; do
    peaktype="binnedPeak"
    prefix=$(basename $binnedPeaks | cut -d '.' -f 1 | sed 's/allChIPCounts_//g')

    prefix="${chip}_${peaktype}_consensusPeaks_${prefix}"

    ## First part
    mergecols=`seq 2 9 | tr '\n' ','`
    expandparam=''



    ## Second part
    # check if the file exists or it was created with a previous peaks version 
    boolAnotMatr=${outpath}/${prefix}.boolean.annotatePeaks.txt
    fileNotExistOrOlder "${boolAnotMatr}" "${binnedPeaks}"
    # this outputs analyse as yes or no in lowercase
    if [[ ${analyse} == "yes" ]]; then

        annotatePeaks.pl \
                ${binnedPeaks} \
                ${speciesGenome} \
                -gid \
                -cpu ${SLURM_CPUS_PER_TASK} \
                -annStats ${outpath}/${prefix}.annotateStats.txt \
                > ${outpath}/${prefix}.annotatePeaks.txt

        cut -f2- ${outpath}/${prefix}.annotatePeaks.txt | \
            awk 'NR==1; NR > 1 {print $0 | "sort -T '.' -k1,1 -k2,2n"}' | \
            cut -f6- > ${outpath}/tmp.txt
        paste ${binnedPeaks} \
            ${outpath}/tmp.txt > ${boolAnotMatr}
    fi
done

echo -e "consensus peak annotations - Finished ---------------------\n"
