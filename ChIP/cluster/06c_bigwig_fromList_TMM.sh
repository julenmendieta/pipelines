#!/bin/bash
# -*- ENCODING: UTF-8 -*-

##===============================================================================
## SLURM VARIABLES
#SBATCH --job-name=bigwig_TMM
#SBATCH --cpus-per-task=12
#SBATCH --mem=5G
#SBATCH --time=05:00:00
#SBATCH -p short
#SBATCH -o /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.out  
#SBATCH -e /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.err 

# HOW TO RUN ME
# add in the list bellow the bam file names (but the .bam termination)
#sbatch /home/jmendietaes/programas/pipelines/ChIP/cluster/06c_bigwig_fromList_TMM.sh


basePath="/home/jmendietaes/data/2021/chip/allProcessed"
# Where to look for bam files and where to store output tree
bamsPath="${basePath}/bamfiles/valid/MEP_vs_GMP"
outpath=${basePath}"/furtherAnalysis/subsampled_noIgG"

# get list of bamfiles to scale
allbams=$(find ${bamsPath}/*bam -printf "%f\n" | \
            tr '\n' ' ')

# path for the location of the pipeline scripts
scriptsPath="/home/jmendietaes/programas/pipelines"

## load modules
export PATH="/home/jmendietaes/programas/miniconda3/envs/Renv/bin:$PATH"
export PATH="/home/jmendietaes/programas/miniconda3/bin:$PATH"


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

############################################
# Scaling factors from allChipCell consensus peaks
############################################

# featureCount files start with this columns before the bam read counts
#Geneid	Chr	Start	End	Strand	Length
# prevFields=6
# lengthCol=5

# featureCpath=${outpath}/peakCalling/MACS2/consensusPeaks/featureCounts
# chip="allmerged"
# for peaktype in narrowPeak broadPeak; do
#     prefix="${chip}_${peaktype}_consensusPeaks"
#     featureOut=${featureCpath}/${prefix}.featureCounts.txt
#     scalingVals=${featureCpath}/${prefix}.featureCounts.TMMscaling.txt

#     Rscript ${scriptsPath}/ChIP/cluster/06_NR_getTMMscaling.R \
#                         -i ${featureOut} \
#                         -o ${scalingVals} \
#                         -l ${lengthCol} \
#                         -d ${prevFields}

# done


############################################
# BigWig normalisation
############################################
echo -e "Starting BigWigs scaling ---------------------------\n"

# Create output dir
if [ ! -e ${outpath}/bigWig_TMM ]; then
    mkdir -p ${outpath}/bigWig_TMM
fi

# for now we only work with narrow peaks
chip="allmerged"
for peaktype in narrowPeak; do
    # Create output dir
    if [ ! -e ${outpath}/bigWig_TMM/${peaktype} ]; then
        mkdir -p ${outpath}/bigWig_TMM/${peaktype}
    fi

    prefix="${chip}_${peaktype}_consensusPeaks"
    scalingVals=${featureCpath}/${prefix}.featureCounts.CPM.TMMscale.txt

    # check content of eleventh line of step control file
    for fi in ${allbams}; do
        echo ${fi}
        
        fi=(${fi//\.bam/ }); fi=${fi[0]}

        # get string of interest from name
        id1=(${fi//\./ }); id1=${id1[0]}
        id1=(${id1//-/ }); id1=${id1[0]}
        id1=(${id1//_/ }); id1="${id1[0]}_${id1[1]}"
        
        scalingFactor=$(grep ${id1} ${scalingVals} | cut -f 2)

        bamPath="${bamsPath}/${fi}.bam"
        bigWigOut2="${outpath}/bigWig_TMM/${peaktype}/${fi}.TMM.bw"

        # check if the file exists of it was created with a previous bam version 
        fileNotExistOrOlder "${bigWigOut2}" "${bamPath}"
        if [[ ${analyse} == "yes" ]]; then
            bamCoverage --binSize 5 --scale ${scalingFactor} \
                -b ${bamPath} -of bigwig \
                -o ${bigWigOut2} --numberOfProcessors $SLURM_CPUS_PER_TASK
        fi

        echo -e "BigWigs scaling - done ---------------------------------------------\n"
    done  
done


echo -e "END --------------------------------------------------"
