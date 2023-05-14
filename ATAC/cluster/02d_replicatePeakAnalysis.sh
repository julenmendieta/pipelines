#!/bin/bash
# -*- ENCODING: UTF-8 -*-

##===============================================================================
## SLURM VARIABLES
#SBATCH --job-name=replyATAC
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=10:00:00
#SBATCH -p short
#SBATCH -o /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.out  
#SBATCH -e /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.err 

# HOW TO RUN ME
#sbatch /home/jmendietaes/programas/PhD/ATAC/cluster/02d_replicatePeakAnalysis.sh \
#/home/jmendietaes/data/2021/ATAC/allProcessed/bamfiles/valid/mergedReplicates/08_paperChip \
#/home/jmendietaes/data/2021/ATAC/allProcessed/furtherAnalysis/08_paperChip 

# NOTES
# As it is now, the script will compare groups with more than one replicates
#   using the whole cell consensus peak table as input for the comparison
#   coordiantes

# path where we have the replicate files
# (inside we have ${replicatesPath}/bamfiles/valid/)
replicatesPath=$1
#replicatesPath="/home/jmendietaes/data/2021/ATAC/allProcessed/bamfiles/valid/mergedReplicates/02_firstATAC"
# Path to the outpath folder where we run 02b_peakAnalysis_I.sh
# and where we have the folder tree with the consensus peaks
inPath=$2
#inPath="/home/jmendietaes/data/2021/ATAC/allProcessed/furtherAnalysis/02_firstATAC"

# If we want to include bamfiles with no replicates in the sae batch
extraBamsPath="/home/jmendietaes/data/2021/ATAC/allProcessed/bamfiles/valid/05-2_laura-NTC"


# State how to merge chip files (appart from whole merge)
# Options are: chip, cellfirst, and cell. That state for
# second name slot in between "_", first name slot in between "_"
# and removing all after "-", and whole first name slot in between "_"
mergeBy="cell"

# Set to TRUE if you want to compare only by batches
#byBatch=TRUE

# Coma separated string with all posible control IDs 
# Used for batch corrected analysis of same KOs
# Set to "no" to not do it
#posibleControls="NTC,WT,NTC0005,NTC5,V12h"
#posibleControls="no"

# extend variables
bamsPath="${replicatesPath}"
outpath=${inPath}"/replicateAnalysis"



# path for the location of the pipeline scripts
scriptsPath="/home/jmendietaes/programas/PhD"

allbams=$(find ${bamsPath}/*bam -printf "${bamsPath}/%f\n" | \
            { grep -v -e "_input" -v -e "_IgG" || :; })

# then extract chip or cell from them
if [[ ${mergeBy} == "chip" ]]; then
    echo "Merging by first element in second slot"
    allChip=$(\
    for filename in ${allbams}; do 
        filename=$(basename ${filename})
        mapLib=(${filename//_/ }); 
        mapLib=${mapLib[1]}; 
        mapLib=(${mapLib//-/ }); 
        echo ${mapLib[0]}; done | \
        grep -v "_input" | grep -v "_IgG" | sort | uniq)

    mergeGroups=${allChip}
elif [[ ${mergeBy} == "cellfirst" ]]; then
    echo "Merging by first element in cell slot"
    allfirstCells=$(\
    for filename in ${allbams}; do 
        filename=$(basename ${filename})
        mapLib=(${filename//_/ }); 
        mapLib=${mapLib[0]}; 
        mapLib=(${mapLib//-/ }); 
        echo ${mapLib[0]}; done | sort | uniq)
    mergeGroups=${allfirstCells}
elif [[ ${mergeBy} == "cell" ]]; then
    echo "Merging by whole cell slot"
    allCells=$(\
    for filename in ${allbams}; do 
        filename=$(basename ${filename})
        mapLib=(${filename//_/ }); 
        mapLib=${mapLib[0]}; 
        echo ${mapLib}; done | sort | uniq)
    mergeGroups=${allCells}
else
    echo "ERROR: Merge-by method not recognised"
    exit 1
fi

#SLURM_CPUS_PER_TASK=6


## load modules
export PATH="/home/jmendietaes/programas/miniconda3/envs/Renv/bin:$PATH"
export PATH="/home/jmendietaes/programas/miniconda3/bin:$PATH"

module load BEDTools/2.27.1-foss-2018b
# fore featureCounts
#module load Subread/1.6.3-foss-2018b

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
                echo $1" older than"${tfile}
            fi
        done
    fi
}

function join_by { local IFS="$1"; shift; echo "$*"; }

###########################################################
# Count reads in consensus peaks with featureCounts
##########################################################
# We will use the generated .saf files in the within chip consensus 
# peaks as the input coordenates

# create new folders
featureCpath=${outpath}/featureCounts
featureCpath_batch=${outpath}/featureCounts_batchCorrect
if [ ! -e ${featureCpath} ]; then
    mkdir -p ${featureCpath}
fi
if [[ ${posibleControls} != "no" ]]; then
    mkdir -p ${featureCpath_batch}
fi


cd ${featureCpath}

echo -e "Starting consensus same-chip featureCounts -----------------------\n"

# get only chips for which we have info of more than one file
chipCheck=$(\
for chip in ${mergeGroups}; do
    #chipFiles=$(echo $allbams | grep -o "\w*${chip}[A-Za-z0-9_\.\-]*")
    chipFiles=$(for a in $allbams; do basename $a | \
                            { grep -e "${chip}-\|${chip}_" || :; } ; done)
    
    cells=$(\
    for filename in ${chipFiles}; do 
        filename=$(basename ${filename})
        mapLib=(${filename//_/ }); 
        echo ${mapLib[0]}; done | sort)

    ncell=$(echo $cells | wc -w)
    if [ "${ncell}" -gt 1 ]; then
        echo $chip
    fi
done)
chipCheck=$(echo $chipCheck | sort)

separator="\|" 
regex="$( printf "${separator}%s" ${chipCheck} )"
regex="${regex:${#separator}}" # remove leading separator
separator="_" 
toCompare="$( printf "${separator}%s" ${chipCheck} )"
toCompare="${toCompare:${#separator}}" # remove leading separator
toCompare=$(echo $toCompare | sed 's/_/-/g')

for peaktype in broadPeak; do
    prefix="${toCompare}_${peaktype}_consensusPeaks"
    coordFiles=$(find -L ${inPath}/peakCalling/MACS2/consensusPeaks/ATAC_${peaktype}_consensusPeaks.saf \
                            -maxdepth 1  -type f ! -size 0 )
    # get a list with the bam files 
    #bamfiles=$(echo $allbams | \
    #                 grep -o "[A-Za-z0-9_\/\.\-]*\w*${chip}[A-Za-z0-9_\.\-]*" | \
    #                 tr '\n' ' ')
    bamfiles=$(echo $allbams | tr ' ' '\n' | \
                        { grep -e ${regex} || :; })

    # check if the file exists or it was created with a previous bam version 
    featureOut=${featureCpath}/${prefix}.featureCounts.txt
    fileNotExistOrOlder "${featureOut}" "${bamfiles}"
    # this outputs analyse as yes or no in lowercase
    if [[ ${analyse} == "yes" ]]; then
        # this only for the consensus between replicates, here no sense
        featureCounts \
                -F SAF \
                -O \
                --fracOverlap 0.2 \
                -T ${SLURM_CPUS_PER_TASK} \
                -p --donotsort \
                -a ${inPath}/peakCalling/MACS2/consensusPeaks/ATAC_${peaktype}_consensusPeaks.saf \
                -o ${featureOut} \
                ${bamfiles} 
    fi
done

echo -e "Consensus same-chip featureCounts - Finished ----------------------\n"


###########################################################
# Differential analysis with DESeq2
###########################################################

if [ ! -e ${outpath}/DESeq2/ ]; then
    mkdir -p ${outpath}/DESeq2/
    mkdir -p ${outpath}/DESeq2_batchCorrect/
fi

for peaktype in broadPeak; do
    prefixF="${toCompare}_${peaktype}_consensusPeaks"

    # first input is featureCounts file
    featureOut=${featureCpath}/${prefixF}.featureCounts.txt

    # check if the file exists or it was created with a previous featureCounts version 
    fileNotExistOrOlder "${outpath}/DESeq2/${toCompare}_${peaktype}/${toCompare}_${peaktype}_DESeq2.log" "${featureOut}"
    if [[ ${analyse} == "yes" ]]; then
        
        
        # Run all by batch
        # Rscript ${scriptsPath}/ATAC/cluster/02_NR_featurecounts_deseq2.r \
        #         --featurecount_file ${featureOut} \
        #         --bam_suffix '.sort.rmdup.rmblackls.rmchr.Tn5.bam' \
        #         --outdir ${outpath}/DESeq2/${toCompare}_${peaktype}/ \
        #         --outprefix $prefix \
        #         --outsuffix '' \
        #         --cores ${SLURM_CPUS_PER_TASK} \
        #         --bybatch ${byBatch}


        # Run merging by ko and with batch corrections
        Rscript ${scriptsPath}/ATAC/cluster/02_NR_featurecounts_deseq2_joinBatches.r \
                --featurecount_file ${featureOut} \
                --bam_suffix '.sort.rmdup.rmblackls.rmchr.Tn5.bam' \
                --outdir ${outpath}/DESeq2_batchCorrect/${toCompare}_${peaktype}/ \
                --outprefix $prefix \
                --outsuffix '' \
                --cores ${SLURM_CPUS_PER_TASK}
    fi
done


# Gather all data
# cd ${outpath}/DESeq2/
# mkdir -p ${outpath}/DESeq2/gatheredDESeq
# for chip in *Peak; do 
#     for compare in ${chip}/*vs*; do 
#         cells=$(basename ${compare}); 
#         cp ${compare}/${cells}.deseq2.results.txt gatheredDESeq/${chip}_${cells}.deseq2.results.txt ; 
#     done;
# done

cd ${outpath}/DESeq2_batchCorrect/
mkdir -p ${outpath}/DESeq2_batchCorrect/gatheredDESeq
for chip in *Peak; do 
    for compare in ${chip}/*vs*; do 
        cells=$(basename ${compare}); 
        cp ${compare}/${cells}.deseq2.results.txt gatheredDESeq/${chip}_${cells}.deseq2.results.txt ; 
    done;
done
