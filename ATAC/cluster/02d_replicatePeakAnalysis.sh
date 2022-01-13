#!/bin/bash
# -*- ENCODING: UTF-8 -*-

##===============================================================================
## SLURM VARIABLES
#SBATCH --job-name=replyPeaks
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH -p short
#SBATCH -o /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.out  
#SBATCH -e /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.err 

# HOW TO RUN ME
#sbatch /home/jmendietaes/programas/PhD/ChIP/cluster/02d_replicatePeakAnalysis.sh \
#/home/jmendietaes/data/2021/chip/allProcessed \
#/home/jmendietaes/data/2021/chip/allProcessed/furtherAnalysis/subsampled_noIgG \


# path where we have the folder structure for chip analysis 
# (inside we have ${basePath}/bamfiles/valid/)
basePath=$1
#basePath="/home/jmendietaes/data/2021/chip/allProcessed"
# Path to the outpath folder where we run 02b_peakAnalysis_I.sh
# and where we have the folder tree with the consensus peaks
inPath=$2
#inPath="${basePath}/furtherAnalysis/subsampled_noIgG"

# path to the header files
extraFilePath="/home/jmendietaes/data/2021/chip/analysisFiles"
# path for the location of the pipeline scripts
scriptsPath="/home/jmendietaes/programas/PhD"

# extend variables
bamsPath="${basePath}/bamfiles/valid/mergedReplicates"
outpath=${basePath}"/furtherAnalysis/replicateAnalysis"


allbams=$(find ${bamsPath}/*bam -printf "${bamsPath}/%f\n" | \
            { grep -v -e "_input" -v -e "_IgG" || :; })
allChip=$(\
    for filename in ${allbams}; do 
        filename=$(basename ${filename})
        mapLib=(${filename//_/ }); 
        mapLib=${mapLib[1]}; 
        mapLib=(${mapLib//-/ }); 
        echo ${mapLib[0]}; done | sort | uniq)
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

###########################################################
# Count reads in consensus peaks with featureCounts
##########################################################
# We will use the generated .saf files in the within chip consensus 
# peaks as the input coordenates

# create new folders
featureCpath=${outpath}/featureCounts
if [ ! -e ${featureCpath} ]; then
	mkdir -p ${featureCpath}
fi
cd ${featureCpath}

echo -e "Starting consensus same-chip featureCounts -----------------------\n"

# get only chips for which we have info of more than one file
chipCheck=$(\
for chip in ${allChip}; do
    chipFiles=$(echo $allbams | grep -o "\w*${chip}[A-Za-z0-9_\.\-]*")

    cells=$(\
    for filename in ${chipFiles}; do 
        filename=$(basename ${filename})
        mapLib=(${filename//_/ }); 
        echo ${mapLib[0]}; done | sort | uniq)

    ncell=$(echo $cells | wc -w)
    if [ "${ncell}" -gt 1 ]; then
        echo $chip
    fi
done)

for chip in ${chipCheck}; do
    for peaktype in narrowPeak broadPeak; do
        prefix="${chip}_${peaktype}_consensusPeaks"
        coordFiles=$(find -L ${inPath}/peakCalling/MACS2/consensusPeaks/bySameChip/${chip}_${peaktype}_consensusPeaks.saf \
                                -maxdepth 1  -type f ! -size 0 )
        # get a list with the bam files 
        bamfiles=$(echo $allbams | grep -o "[A-Za-z0-9_\/\.\-]*\w*${chip}[A-Za-z0-9_\.\-]*" | tr '\n' ' ')

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
                    -a ${inPath}/peakCalling/MACS2/consensusPeaks/bySameChip/${prefix}.saf \
                    -o ${featureOut} \
                    ${bamfiles} 
        fi
    done
done

echo -e "Consensus same-chip featureCounts - Finished ----------------------\n"


###########################################################
# Differential analysis with DESeq2
###########################################################

if [ ! -e ${outpath}/DESeq2/ ]; then
	mkdir -p ${outpath}/DESeq2/
fi

for chip in ${chipCheck}; do
    for peaktype in narrowPeak broadPeak; do
        prefixF="${chip}_${peaktype}_consensusPeaks"

        # first input is featureCounts file
        featureOut=${featureCpath}/${prefixF}.featureCounts.txt

        # get compared cells to add them at name
        cells=$(head -n 2 ${featureOut} | tail -n 1 \
                | grep -o "[\/]\w*_*${chip}")
        cells=`for cell in ${cells}; do 
            mapLib=(${cell//_/ }); 
            mapLib=${mapLib[0]}; 
            echo ${mapLib} | sed 's/\///g'; done | sort | uniq | \
            tr '\n' '-'`
        cells=${cells::-1}
        prefix="${chip}_${cells}_${peaktype}_DESeq2" 

        Rscript ${scriptsPath}/ChIP/cluster/02_NR_featurecounts_deseq2.r \
                --featurecount_file ${featureOut} \
                --bam_suffix '.sort.rmdup.rmblackls.rmchr.bam' \
                --outdir ${outpath}/DESeq2/${chip}_${peaktype}/ \
                --outprefix $prefix \
                --outsuffix '' \
                --cores ${SLURM_CPUS_PER_TASK}
    done
done

