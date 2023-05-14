#!/bin/bash
# -*- ENCODING: UTF-8 -*-

##===============================================================================
## SLURM VARIABLES
#SBATCH --job-name=mergeSeqRep
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=00-24:00:00
#SBATCH -p short
#SBATCH -o /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.out  
#SBATCH -e /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.err 
##SBATCH --dependency=afterany:793158


##SBATCH --mail-type=END
##SBATCH --mail-user=user@mail.es
# HOW TO RUN ME
# sbatch /home/jmendietaes/programas/PhD/ChIP/cluster/11a_mergeSeqRep.sh 


# PURPOSE
# To merge sequencing replicates bam files and run duplicate removal again

## path where the sequencing replicates are stored
inBam="/home/jmendietaes/data/2021/chip/allProcessed/bamfiles/valid/mergedSeqRep"
# The file with the bams to merge is in this folder as mergedAndSubs.txt
# lines not starting with # should have all the bams to merge separated
# by a space, with last name being the output merged bam name

## path where the output merged bam files will be stored
outBam="/home/jmendietaes/data/2021/chip/allProcessed/bamfiles/valid"
## path where the output merged Bigwig files will be stored
outBw="/home/jmendietaes/data/2021/chip/allProcessed/BigWig/valid"

# cells to be counted at the time to merge (separated by \|)
#subsampleCells="DM_\|Mye_"

# Set to "yes" if you want to remove the temporal bam files created at the time
# to remove duplicates
removeTemp="yes"

mergeFile="${inBam}/mergedAndSubs.txt"
##===============================================================================
wigToBigWig="/home/jmendietaes/programas/binPath/wigToBigWig"
picardPath='/home/jmendietaes/programas/picard/picard.jar'
## Required Software
module load SAMtools/1.12-GCC-10.2.0
module load Java/1.8.0_192
#module load picard/2.18.17-Java-1.8

##===============================================================================

# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

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
nCPU=${SLURM_CPUS_PER_TASK}

# Create temp output dir
if [ ! -e ${inBam}/temp ]; then
    mkdir -p ${inBam}/temp
fi


cd ${inBam}
#checkBams=$(ls *bam | grep "${subsampleCells}")
mergeLines=$(grep -v "#" ${mergeFile} | tr ' ' ',')

for merge in ${mergeLines}; do
    echo "Merging new files"
    echo $merge
    ##############################
    # File merging
    ##############################
    echo -e "Starting File merging ------------------------------------------------ \n"

    # Get files to merge and output name separately
    merge=(${merge//\,/ }); 
    mergeB=$(echo "${merge[@]::${#merge[@]}-1}")
    mergeOut=${merge[-1]}
    mergeOutDup="withDup_${mergeOut}"

    # check if the final bam exists or it was created with a previous bam version
    # This in the main output directory 
    fileNotExistOrOlder "${outBam}/${mergeOut}" "${mergeB}"
    # this outputs analyse as yes or no in lowercase
    if [[ ${analyse} == "yes" ]]; then
        # check if the final bam exists or it was created with a previous bam version
        # This in the merged Replicates directory 
        fileNotExistOrOlder "${outBam}/mergedReplicates/${mergeOut}" "${mergeB}"
        if [[ ${analyse} == "yes" ]]; then

            # check if the merged un-dedup bam exists or it was created with a previous bam version 
            fileNotExistOrOlder "${inBam}/temp/${mergeOutDup}" "${mergeB}"
            # this outputs analyse as yes or no in lowercase
            if [[ ${analyse} == "yes" ]]; then

                # merge bam files
                samtools merge --threads ${nCPU} ${inBam}/temp/${mergeOutDup} ${mergeB}
                #mv ${outName}* ../
            else
                echo -e "File merging - already done before  ---------------------------- \n"
            fi
        else
            echo -e "File merging - already done before  ---------------------------- \n"
        fi

    else
        echo -e "File merging - already done before  ---------------------------- \n"
    fi
    ##############################
    # Dupicates marking and removal
    ##############################
    echo -e "Starting mark and remove duplicates ------------------------------------------------ \n"
    # Samtools flags to remove duplicates + unmapped reads:
    #   a) -F 1028 => read unmapped // duplicate
    #   b) -F 1804 => read unmapped // mate unmapped // not primary alignment // read failing platform // duplicate

    mergeOutMarkDup="markDup_${mergeOut}"
    
    # check if the final bam exists or it was created with a previous bam version
    # This in the main output directory 
    fileNotExistOrOlder "${outBam}/${mergeOut}" "${mergeB}"
    # this outputs analyse as yes or no in lowercase
    if [[ ${analyse} == "yes" ]]; then

        # check if the final bam exists or it was created with a previous bam version
        # This in the merged Replicates directory 
        fileNotExistOrOlder "${outBam}/mergedReplicates/${mergeOut}" "${mergeB}"
        if [[ ${analyse} == "yes" ]]; then

            # Identify duplicated reads by flaging them in a new ioutput BAM file
            java -Xmx24G -jar ${picardPath} MarkDuplicates \
                I=${inBam}/temp/${mergeOutDup} O=${inBam}/temp/${mergeOutMarkDup} \
                METRICS_FILE=${inBam}/temp/${mergeOutMarkDup}.txt
            # remove reads marked as duplicated
            samtools view -o ${outBam}/${mergeOut} -@ ${nCPU} -bh \
                -F 1804 ${inBam}/temp/${mergeOutMarkDup}
            samtools index ${outBam}/${mergeOut}

            if [[ $removeTemp == 'yes' ]] ; then
                rm ${inBam}/temp/${mergeOutDup}
                rm ${inBam}/temp/${mergeOutMarkDup}
                rm ${inBam}/temp/${mergeOutMarkDup}.txt
            fi
            echo -e "Remove duplicates - done ------------------------------------------------ \n"
            bamForBw="${outBam}/${mergeOut}"
        else
            echo -e "Remove duplicates - already done and merged before ---------- \n"
            bamForBw="${outBam}/mergedReplicates/${mergeOut}"
        fi

    else
        echo -e "Remove duplicates - already done before ------------------------- \n"
        bamForBw="${outBam}/${mergeOut}"
    fi 

    ###################################
    # BigWigs 
    #################################
    # deeptools bamCoverage
    # Normalize by CPM (This is the scaled bigWig we will use)
    echo -e "Starting BigWigs --------------------------------------------------\n"
    bigWigOut=(${mergeOut//bam/ }); 
    bigWigOut=${bigWigOut[0]}
    bigWigOut="${outBw}/${bigWigOut[0]}norm.bw"

    # check if the output exists of it was created with a previous bam version 
    fileNotExistOrOlder "${bigWigOut}" "${bamForBw}"
    # this outputs analyse as yes or no in lowercase
    if [[ ${analyse} == "yes" ]]; then

        # check if the final bam exists or it was created with a previous bam version
        # This in the merged Replicates directory 
        fileNotExistOrOlder "${outBw}/mergedReplicates/${bigWigOut[0]}norm.bw" "${bamForBw}"
        if [[ ${analyse} == "yes" ]]; then

            bamCoverage --binSize 5 --normalizeUsing CPM --exactScaling \
                -b ${bamForBw} -of bigwig \
                -o ${bigWigOut} --numberOfProcessors ${nCPU}
            
            echo -e "BigWig norm - done ---------------------------------------------\n"
        else
            echo -e "BigWig norm - already done before ------------------------------\n"
        fi
    else
        echo -e "BigWig norm - already done before ------------------------------\n"
    fi


done