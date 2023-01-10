#!/bin/bash
# -*- ENCODING: UTF-8 -*-

##===============================================================================
# HOW TO RUN ME
#bash /home/julen/programas/PhD/ChIP/csaw/02b_specificConsensus.sh 

# OBJECTIVE
# Annotate peakfiles in a folder

basePath="/scratch/julen/ChIP/allData/04_subsamplingNoIgG/outdata/csaw_ainhoa"
bamsPath="/scratch/julen/ChIP/bamFiles/paperAinhoa"
binnedCoord="/scratch/julen/ChIP/allData/04_subsamplingNoIgG/outdata/csaw_ainhoa/binnedPeaks/allChIPCounts_DM-Mye-GMP-Mono_IgG-flt.tsv"
outpath="${basePath}/extraConsensus"
scriptsPath="/home/julen/programas/PhD"

nCPU=16

# species shortcut for MACS
speciesGenome="mm10"

## load modules
export PATH="/home/julen/miniconda3/envs/python3/bin:$PATH"
export PATH="/home/julen/miniconda3/envs/onlyR/bin:$PATH"


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


###########################################################
# Count reads in consensus peaks with featureCounts
##########################################################

featureCpath="${outpath}/featureCounts"
if [ ! -e ${featureCpath} ]; then
    mkdir -p ${featureCpath}
fi
cd ${featureCpath}

echo -e "Starting consensus featureCounts -----------------------\n"



prefix="fromBinConsensus"

# get list of bamfiles and labels
bamfiles=$(find -L ${bamsPath}/*bam -printf "${bamsPath}/%f " \
            | tr '\n' ' ')
fileLabels=$(for f in $bamfiles; do echo ${f##*/} |sed "s/.sort.rmdup.rmblackls.rmchr.bam//g"; done)

# we create the saff file
echo -e "GeneID\tChr\tStart\tEnd\tStrand" > ${featureCpath}/consensus.saf
awk '{print $4"\t"$1"\t"$2"\t"$3"\t+"}' $binnedCoord | tail -n +2 >> ${featureCpath}/consensus.saf

# check if the file exists or it was created with a previous bams list 
featureOut=${featureCpath}/${prefix}.featureCounts.txt
fileNotExistOrOlder "${featureOut}" "${bamsPath}/*bam"
# this outputs analyse as yes or no in lowercase
if [[ ${analyse} == "yes" ]]; then
    # this only for the consensus between replicates, here no sense
    featureCounts \
            -F SAF \
            -O \
            --fracOverlap 0.2 \
            -T ${nCPU} \
            -p --donotsort \
            -a ${featureCpath}/consensus.saf \
            -o ${featureOut} \
            ${bamfiles}
fi
echo -e "Consensus featureCounts - Finished ----------------------\n"



###########################################################
# Read count to CPM
###########################################################

echo -e "Starting consensus CPM -----------------------\n"

# featureCount files start with this columns before the bam read counts
#Geneid	Chr	Start	End	Strand	Length
prevFields=6
lengthCol=5

# set the output file path and copy the content of original
featureCPM=${featureCpath}/${prefix}.featureCounts.CPM.txt

echo -e "FileName\tSampleName\tCellType\tStatus" > ${featureCpath}/sampleInfo.txt
for bam in ${bamfiles}; do
    sname=`basename $bam  | sed 's/.sort.rmdup.rmblackls.rmchr.bam//g'`
    #sname=$bam
    cell=(${sname//_/ }) ; cell=${cell[0]}
    status=(${sname//_/ }) ; status=${status[1]}; status=(${status//-/ }); status=${status[0]}
    echo -e ${bam}"\t"${sname}"\t"${cell}"\t"${status} >> ${featureCpath}/sampleInfo.txt
done

# check if the file exists or it was created with a previous featureCounts version 
fileNotExistOrOlder "${featureCPM}" "${featureCpath}/${prefix}.featureCounts.txt"
# this outputs analyse as yes or no in lowercase
if [[ ${analyse} == "yes" ]]; then
    # this only for the consensus between replicates, here no sense
    Rscript ${scriptsPath}/ChIP/cluster/02_NR_metric_CPMfromFeatureCounts.r \
                -i ${featureCpath}/${prefix}.featureCounts.txt \
                -o ${featureCPM} \
                -s ${featureCpath}/sampleInfo.txt \
                -l ${lengthCol} \
                -d ${prevFields} 
                    

fi

echo -e "Consensus CPM - Finished ----------------------\n"


