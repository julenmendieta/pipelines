#!/bin/bash
# -*- ENCODING: UTF-8 -*-

##===============================================================================
## SLURM VARIABLES
#SBATCH --job-name=SSP
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=10:00:00
#SBATCH -p short
#SBATCH -o /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.out  
#SBATCH -e /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.err 


# HOW TO RUN ME
#sbatch /home/jmendietaes/programas/pipelines/ChIP/cluster/08a_SSP.sh \
#/home/jmendietaes/data/2021/chip/allProcessed/bamfiles/valid \
#/home/jmendietaes/data/2021/chip/allProcessed/SSP \
#/home/jmendietaes/referenceGenomes/mm10_reordered/mm10.reordered.sizes.shortTab

# OBJECTIVE
# Run SSP to get metrics of ChIP-seq validity (enough sequencing etc)
# Sensitive and robust assessment of ChIP-seq read distribution using 
#a strand-shift profile Ryuichiro Nakato* and Katsuhiko Shirahige

# Path to bam files (we assume non-merged replicates are in that path 
#inside mergedReplicates folder)
bamsPath=$1

# Path were to store output
outpath=$2

# Path to the tab separated chromosome sizes file (only the non filtered ones)
genomeTable=$3

# path to SSP singularity path
sspPath="/home/jmendietaes/programas/SSP/ssp_drompa.img"

## load modules
export PATH="/home/jmendietaes/programas/miniconda3/envs/Renv/bin:$PATH"
module load Boost/1.77.0-GCC-11.2.0
module load Singularity/3.10.2

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

#############################
# Run SSP
#############################

# First we need the paths to all bam files
allbams=$(find ${bamsPath}/*bam -printf "${bamsPath}/%f\n" | \
            tr '\n' ' ')

allbams2=$(find ${bamsPath}/mergedReplicates/*bam -printf "${bamsPath}/mergedReplicates/%f\n" | \
            tr '\n' ' ')
allbams=$(echo $allbams $allbams2)

for bam in ${allbams}; do
    id=$(basename $bam); id=(${id//./ }); id=${id[0]}
    outfile=${id}.stats.txt

    # check if the file exists of it was created with a previous bam version 
    fileNotExistOrOlder "${outfile}" "${bam}"
    # this outputs analyse as yes or no in lowercase

    if [[ ${analyse} == "yes" ]]; then

        singularity exec ${sspPath} ssp -i ${bam} -o ${id} --odir ${outpath} \
                                        --pair -f BAM --gt ${genomeTable} \
                                        -p ${SLURM_CPUS_PER_TASK}
    fi
done

statsFiles=`find ${outpath}/*stats.txt`
statsFile1=(${statsFiles// / }); 
statsFile1=${statsFile1[0]}
head -n 1 $statsFile1 > ${outpath}/allStats_SSP.txt
for fi in ${statsFiles}; do
     tail -n 1 ${fi} >> ${outpath}/allStats_SSP.txt
done
