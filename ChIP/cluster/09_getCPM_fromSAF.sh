#!/bin/bash
# -*- ENCODING: UTF-8 -*-

##===============================================================================
## SLURM VARIABLES
#SBATCH --job-name=cpmInFocusCoord
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --time=10:00:00
#SBATCH -p short
#SBATCH -o /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.out  
#SBATCH -e /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.err 


# HOW TO RUN ME
#sbatch /home/jmendietaes/programas/pipelines/ChIP/cluster/09_getCPM_fromSAF.sh 

# OBJECTIVE
# get CPM values from bamfiles in path at consensus saf coordinta file

# For my paper with ATAC
peaktype="broadPeak"
#peaktype="narrowPeak"
# Path to coordinate file
safIn="/home/jmendietaes/data/2021/chip/allProcessed/furtherAnalysis/08_projectRestart/csaw/DM_RPKM.saf"
#safIn="/home/jmendietaes/data/2021/chip/allProcessed/furtherAnalysis/08_projectRestart/csaw/MACSconsensus_DM-Mye_WS-200_SPC-200_IgG-flt.saf"
#safIn="/home/jmendietaes/data/2021/chip/allProcessed/furtherAnalysis/15_ATACinTF/specificConsensus/DM_TFcentred.saf"
#safIn="/home/jmendietaes/data/2021/chip/allProcessed/furtherAnalysis/15_ATACinTF/specificConsensus/CFaroundATAC.saf"
#safIn="/home/jmendietaes/data/2021/ATAC/allProcessed/furtherAnalysis/08c_inUMAP/peakCalling/MACS2/consensus/allmerged_broadPeak_consensusPeaks.saf"
#safIn="/home/jmendietaes/data/2021/ATAC/allProcessed/furtherAnalysis/08d_chipOnATAC/peakCalling/MACS2/consensusPeaks/DM_allpeaks.saf"
#safIn="/home/jmendietaes/data/2021/ATAC/allProcessed/furtherAnalysis/08d_chipOnATAC/peakCalling/MACS2/consensusPeaks/DM_allpeaks_4kbext.saf"
# Path where the bams of interest are stored
bamsIn="/home/jmendietaes/data/2021/ATAC/allProcessed/bamfiles/valid/08_paperChip"
#bamsIn="/home/jmendietaes/data/2021/chip/allProcessed/bamfiles/valid/15_CFinTF/"
#bamsIn="/home/jmendietaes/data/2021/chip/allProcessed/bamfiles/valid/16_TFcentred"
#bamsIn="/home/jmendietaes/data/2021/chip/allProcessed/bamfiles/valid/16b_TFcentred_merged"
#bamsIn="/home/jmendietaes/data/2021/ATAC/allProcessed/bamfiles/valid/08c_inUMAP"
#bamsIn="/home/jmendietaes/data/2021/chip/allProcessed/bamfiles/valid/08d_chipOnATAC"
# Path where to output files
#outPath="/home/jmendietaes/data/2021/chip/allProcessed/furtherAnalysis/08_projectRestart/csaw/ATACcounts"
outPath="/home/jmendietaes/data/2021/chip/allProcessed/furtherAnalysis/08_projectRestart/csaw/ATACcountsChunk"
#outPath="/home/jmendietaes/data/2021/chip/allProcessed/furtherAnalysis/15_ATACinTF/specificConsensus"
#outPath="/home/jmendietaes/data/2021/ATAC/allProcessed/furtherAnalysis/08c_inUMAP/counts"
#outPath="/home/jmendietaes/data/2021/chip/allProcessed/furtherAnalysis/08d_chipOnATAC/counts"
# Files first label
#label1="ATAC"
label1="DM_ATAC"
#label1="CFaroundATAC"
#label1="DM_TFcentred"
#label1="DM_TFcentred_merged"
#label1="chipOnATAC"
#label1="chipOnATAC_4kbext"
# Bam file end string to be removed from name
#bamEnd=".sort.rmdup.rmblackls.rmchr.Tn5.bam"
bamEnd=".sort.rmdup.*"
#bamEnd=".sort.rmdup.rmblackls.rmchr.bam"

# For Ainhoa paper
# safIn="/home/jmendietaes/data/2021/chip/allProcessed/furtherAnalysis/07_leukPaper/peakCalling/MACS2/consensusPeaks/allmerged_narrowPeak_consensusPeaks.saf"
# peaktype="narrowPeak"
# bamsIn="/home/jmendietaes/data/2021/chip/allProcessed/bamfiles/valid/07_leukPaper3_addings"
# outPath="/home/jmendietaes/data/2021/chip/allProcessed/furtherAnalysis/07_leukPaper/addMany"
# label1="leukAddings"

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
# the second input file/s
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

cd ${outPath}

echo -e "Starting consensus featureCounts -----------------------\n"


prefix="${label1}_${peaktype}_consensusPeaks"
bamfiles=$(find -L ${bamsIn}/*.bam -maxdepth 1  -type f ! -size 0)
fileLabels=$(for f in $bamfiles; do echo ${f##*/} |sed "s/${bamEnd}//g"; done)


# check if the file exists or it was created with a previous peaks version 
featureOut=${outPath}/${prefix}.featureCounts.txt
fileNotExistOrOlder "${featureOut}" "${bamfiles} ${safIn}"
# this outputs analyse as yes or no in lowercase
if [[ ${analyse} == "yes" ]]; then
    # this only for the consensus between replicates, here no sense
    featureCounts \
            -F SAF \
            -O \
            --fracOverlap 0.2 \
            -T ${SLURM_CPUS_PER_TASK} \
            -p --donotsort \
            -a ${safIn} \
            -o ${featureOut} \
            ${bamfiles}
fi
echo -e "Consensus featureCounts - Finished ----------------------\n"



###########################################################
# Read count to CPM
###########################################################

echo -e "Starting consensus CPM -----------------------\n"

nbamfiles=$(echo $bamfiles | wc -w)
# featureCount files start with this columns before the bam read counts
#Geneid	Chr	Start	End	Strand	Length
prevFields=6
lengthCol=5

# set the output file path and copy the content of original
featureCPM=${outPath}/${prefix}.featureCounts.CPM.txt

echo -e "FileName\tSampleName\tCellType\tStatus" > ${outPath}/sampleInfo.txt
for bam in ${bamfiles}; do
    sname=`basename $bam  | sed "s/${bamEnd}//g"`
    #sname=$bam
    cell=(${sname//_/ }) ; cell=${cell[0]}
    status=(${sname//_/ }) ; status=${status[1]}; status=(${status//-/ }); status=${status[0]}
    echo -e ${bam}"\t"${sname}"\t"${cell}"\t"${status} >> ${outPath}/sampleInfo.txt
done

# check if the file exists or it was created with a previous featureCounts version 
fileNotExistOrOlder "${featureCPM}" "${outPath}/${prefix}.featureCounts.txt"
# this outputs analyse as yes or no in lowercase
if [[ ${analyse} == "yes" ]]; then
    # this only for the consensus between replicates, here no sense
    Rscript ${scriptsPath}/ChIP/cluster/02_NR_metric_CPMfromFeatureCounts.r \
                -i ${outPath}/${prefix}.featureCounts.txt \
                -o ${featureCPM} \
                -s ${outPath}/sampleInfo.txt \
                -l ${lengthCol} \
                -d ${prevFields} 
                    

fi

echo -e "Consensus CPM - Finished ----------------------\n"

############################################
# Merge CPM with rcounts
############################################
echo -e "Starting info merge ---------------------------\n"

outfile=${outPath}/${prefix}.rcount.CPM.txt


cat ${featureCPM} > ${outPath}/tmp2.txt
headerNames=$(head -n 1 ${outPath}/tmp2.txt)
for h in ${headerNames}; do
    if [[ $h != "interval_id" ]]; then
        sed -i "s/${h}/${h}.cpm/g" ${outPath}/tmp2.txt
    fi
done

# Both files are in the same order
tail -n +2 ${outPath}/${prefix}.featureCounts.txt | \
            sed "s/${bamEnd}/.rcount/g" | \
            sed "s/$(echo $bamsIn/ | sed 's_/_\\/_g')//g" | \
            paste - \
            ${outPath}/tmp2.txt > ${outfile}

# Remove extra path from header
#sed -i "s/$(echo $bamsPath/ | sed 's_/__g')//g" \
#    ${outfile}


echo -e "END --------------------------------------------------"
