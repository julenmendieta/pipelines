#!/bin/bash
# -*- ENCODING: UTF-8 -*-

##===============================================================================
## SLURM VARIABLES
#SBATCH --job-name=prePeak_QC
#SBATCH --cpus-per-task=16
#SBATCH --mem=50G
#SBATCH --time=1-00:00:00
#SBATCH -p medium
#SBATCH -o /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.out  
#SBATCH -e /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.err 


# HOW TO RUN ME
# create file with ChIP names in bams folder
#sbatch /home/jmendietaes/programas/pipelines/ChIP/cluster/02a_prePeak_QC.sh \
#/home/jmendietaes/data/2021/chip/allProcessed \
#/home/jmendietaes/data/2021/chip/analysisFiles

# path where we have the folder structure for chip analysis 
# (inside we have ${basePath}/bamfiles/valid/)
basePath=$1
#basePath="/home/jmendietaes/data/2021/chip/allProcessed"
# Location of the analysis header files downloaded from 
# https://github.com/nf-core/chipseq/tree/master/assets/multiqc
extraFilePath=$2
#extraFilePath="/home/jmendietaes/data/2021/chip/analysisFiles"


# path for the location of the pipeline scripts
scriptsPath="/home/jmendietaes/programas/PhD"

# extend variables
bamsPath="${basePath}/bamfiles/valid"
outpath=${basePath}"/furtherAnalysis"
# Number of genomic bins to use when calculating fingerprint plot (Default: 500000)
fingerprint_bins=500000
allbams=$(find ${bamsPath}/*bam -printf "${bamsPath}/%f ")
allLabels=`for i in $allbams; do basename ${i} | cut -d '.' -f 1 | cut -d '_' -f 1,2,3; done | tr '\n' ' '`
#SLURM_CPUS_PER_TASK=6


## load modules
export PATH="/home/jmendietaes/programas/miniconda3/envs/Renv/bin:$PATH"
export PATH="/home/jmendietaes/programas/miniconda3/bin:$PATH"

module load SAMtools/1.12-GCC-10.2.0

##===============================================================================
# enable debuging displaying line number and content
# export PS4='$LINENO+ '
# set -x 

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
            fi
        done
    fi
}
# #########################################
# # Correlate signal at gene positions
# ########################################

# bed="/home/jmendietaes/data/2021/chip/analysisFiles/hg38Tables_genePos.txt"

# if [ ! -e ${outpath}/chipCorrelation/ ]; then
#     mkdir -p ${outpath}/chipCorrelation/
# fi

# # get all bigwig file names
# allBigW=$(find ${basePath}/BigWig/valid/*bw -printf "${basePath}/BigWig/valid/%f ")

# # Correlate only signal in gene regions 
# computeMatrix scale-regions \
#     --regionsFileName $bed \
#     --scoreFileName $allBigW \
#     --outFileName ${outpath}/chipCorrelation/allChIP.computeMatrix.mat.gz \
#     --outFileNameMatrix ${outpath}/chipCorrelation/allChIP.computeMatrix.vals.mat.tab \
#     --regionBodyLength 1000 \
#     --beforeRegionStartLength 3000 \
#     --afterRegionStartLength 3000 \
#     --skipZeros \
#     --smartLabels \
#     --numberOfProcessors ${SLURM_CPUS_PER_TASK}

# plotProfile --matrixFile ${outpath}/chipCorrelation/allChIP.computeMatrix.mat.gz \
#     --outFileName ${outpath}/chipCorrelation/allChIP.plotProfile.pdf \
#     --outFileNameData ${outpath}/chipCorrelation/allChIP.plotProfile.tab
# plotHeatmap --matrixFile ${outpath}/chipCorrelation/allChIP.computeMatrix.mat.gz \
#     --outFileName ${outpath}/chipCorrelation/allChIP.plotHeatmap.pdf \
#     --outFileNameMatrix ${outpath}/chipCorrelation/allChIP.plotHeatmap.mat.tab



#############################
# Phantompeakqualtools
#############################
# RSC values larger than 1 for all ChIP samples suggest good signal-to-noise
# Low RSC values can be due to failed or poor quality ChIP, low read 
# sequence quality and hence lots of mismappings, shallow sequencing depth 
# or a combination of these. Also, datasets with few binding sites (< 200) 
# help: https://hbctraining.github.io/Intro-to-ChIPseq/lessons/06_combine_chipQC_and_metrics.html
# The red vertical line shows the dominant peak at the true peak shift, 
# with a small bump at the blue vertical line representing the read length
# More than one red line implies weak signal
# https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionV/lessons/CC_metrics_extra.html


# Phantompeakqualtools
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

# BiocManager::install("Rsamtools")
# install.packages("spp", dependencies=TRUE)
#wget https://raw.githubusercontent.com/nf-core/chipseq/master/assets/multiqc/spp_rsc_header.txt
#wget https://raw.githubusercontent.com/nf-core/chipseq/master/assets/multiqc/spp_nsc_header.txt
#wget https://raw.githubusercontent.com/nf-core/chipseq/master/assets/multiqc/spp_correlation_header.txt
#wget https://raw.githubusercontent.com/crazyhottommy/phantompeakqualtools/master/run_spp.R

spp_correlation_header="${extraFilePath}/spp_correlation_header.txt"
spp_nsc_header="${extraFilePath}/spp_nsc_header.txt"
spp_rsc_header="${extraFilePath}/spp_rsc_header.txt"

# create common nsc and rsc file
if [ ! -e ${outpath}/phantompeak/allChIP_spp_nsc.tsv ]; then
    cat $spp_nsc_header > ${outpath}/phantompeak/allChIP_spp_nsc.tsv
    cat $spp_rsc_header > ${outpath}/phantompeak/allChIP_spp_rsc.tsv
fi

# Create output dir
if [ ! -e ${outpath}/phantompeak/extraData ]; then
    mkdir -p ${outpath}/phantompeak/extraData
fi

RUN_SPP="${scriptsPath}/ChIP/cluster/02_NR_run_spp.R"

for bam in $allbams; do
    label=`basename ${bam} | cut -d '.' -f 1 | cut -d '_' -f 1,2,3`

    # check if the file exists of it was created with a previous bam version 
    fileNotExistOrOlder "${outpath}/phantompeak/extraData/${label}.spp.out" "${bam}"
    # this outputs analyse as yes or no in lowercase

    if [[ ${analyse} == "yes" ]]; then

        Rscript -e "library(caTools); source(\"$RUN_SPP\")" -c="${bam}" \
        -savp="${outpath}/phantompeak/${label}.spp.pdf" \
        -savd="${outpath}/phantompeak/${label}.spp.Rdata" \
        -out="${outpath}/phantompeak/${label}.spp.out" -p=${SLURM_CPUS_PER_TASK}

        cp $spp_correlation_header ${outpath}/phantompeak/extraData/${label}_spp_correlation_mqc.tsv

        Rscript -e "load('${outpath}/phantompeak/${label}.spp.Rdata'); \
        write.table(crosscorr\$cross.correlation, \
        file=\"${outpath}/phantompeak/extraData/${label}_spp_correlation_mqc.tsv\", sep=",", \
        quote=FALSE, row.names=FALSE, col.names=FALSE,append=TRUE)"
        rm ${outpath}/phantompeak/${label}.spp.Rdata
        # I will merge this files after the full analysis
        awk -v name=$label -v OFS='\t' '{print name, $9}' ${outpath}/phantompeak/${label}.spp.out >> ${outpath}/phantompeak/allChIP_spp_nsc.tsv
        awk -v name=$label -v OFS='\t' '{print name, $10}' ${outpath}/phantompeak/${label}.spp.out >> ${outpath}/phantompeak/allChIP_spp_rsc.tsv

        # move the files we most surely wont use
        cat ${outpath}/phantompeak/${label}.spp.out >> ${outpath}/phantompeak/extraData/allChIP.spp.out
        mv ${outpath}/phantompeak/${label}.spp.out ${outpath}/phantompeak/extraData/${label}.spp.out
    fi

done



#############################
#FingerPrint
#############################
# did the antibody-treatment enrich sufficiently so that the ChIP 
# signal can be separated from the background signal?
# This tool samples indexed BAM files and plots a profile of cumulative 
# read coverages for each. All reads overlapping a window (bin) of the 
# specified length are counted; these counts are sorted and the cumulative 
# sum is finally plotted.

# Create output dir
if [ ! -e ${outpath}/fingerPrint ]; then
    mkdir -p ${outpath}/fingerPrint
fi

# we will do it per chipped protein and cell against all controls
## First we get all chipped proteins IDs
allChip=`for filename in ${allLabels}; do 
            mapLib=(${filename//_/ }); mapLib=${mapLib[1]}; mapLib=(${mapLib//-/ }); 
            echo ${mapLib[0]}; done | grep -v input | grep -v IgG | sort | uniq`

allCell=`for filename in ${allLabels}; do 
            mapLib=(${filename//_/ }); 
            echo ${mapLib[0]}; done | sort | uniq`

## then we go for each cell
for cell in `echo ${allCell} | tr ' ' '\n' `; do
    echo $cell
    cellChip=`echo ${allLabels} | tr ' ' '\n' | grep ${cell}`
    cellChip=`for filename in ${cellChip}; do 
            mapLib=(${filename//_/ }); mapLib=${mapLib[1]}; mapLib=(${mapLib//-/ }); 
            echo ${mapLib[0]}; done | grep -v input | grep -v IgG | sort | uniq`
    # Dont know why but in submited jobs this might fail when there are no controls
    # so i need a back up plan
    cellControls=`echo ${allbams} | tr ' ' '\n' | grep "${cell}_" | \
                    { grep -e "input" -e "IgG" || :; }`

    # proceed only if we have controls
    if [[ ! $cellControls == "" ]] ; then
        ## Compare all the bamfiles from same cell and chip
        for chip in ${cellChip}; do
            echo $chip
            chipbams=`echo ${allbams} | tr ' ' '\n' | grep "${cell}_${chip}"`
            labels=$(for file in ${chipbams} ${cellControls}; do basename $file | cut -d '.' -f 1 | cut -d '_' -f 1,2,3; done)
            
            # CHECK: might want to change this in the future        
            # for now we use the first control to normalise signal (by sort should be IgG)
            useControl=`echo ${cellControls} | cut -d ' ' -f 1`
            
            outRaw=${outpath}/fingerPrint/${cell}_${chip}_${fingerprint_bins}bp_fingerprint.raw.txt
            # check if the file exists of it was created with a previous bam version 
            fileNotExistOrOlder "${outRaw}" "${chipbams} ${cellControls}"
            # this outputs analyse as yes or no in lowercase

            if [[ ${analyse} == "yes" ]]; then

                plotFingerprint \
                    --bamfiles ${chipbams} ${cellControls} \
                    --plotFile ${outpath}/fingerPrint/${cell}_${chip}_${fingerprint_bins}bp_fingerprint.pdf \
                    --labels ${labels} \
                    --outRawCounts ${outRaw} \
                    --outQualityMetrics ${outpath}/fingerPrint/${cell}_${chip}_${fingerprint_bins}bp_fingerprint.qcmetrics.txt \
                    --skipZeros \
                    --JSDsample ${useControl} \
                    --numberOfProcessors ${SLURM_CPUS_PER_TASK} \
                    --numberOfSamples $fingerprint_bins
            fi
        done
    fi
done


#############################
# Correlations
#############################

if [ ! -e ${outpath}/chipCorrelation/allVSall ]; then
    mkdir -p ${outpath}/chipCorrelation/allVSall
fi


# generate correlation matrix
allBamFiles=`ls  ${bamsPath}/*bam`
labels=$(for file in ${allBamFiles}; do basename $file | cut -d '.' -f 1 | cut -d '_' -f 1,2,3; done)

outMatrix=${outpath}/chipCorrelation/allVSall/chipCorrelations.npz
# check if the file exists of it was created with a previous bam version 
fileNotExistOrOlder "${outMatrix}" "${allBamFiles}"
# this outputs analyse as yes or no in lowercase

if [[ ${analyse} == "yes" ]]; then
    multiBamSummary bins --bamfiles ${allBamFiles} \
                    -o ${outMatrix} \
                    --labels $labels --binSize 10000 \
                    --numberOfProcessors ${SLURM_CPUS_PER_TASK}

    plotCorrelation --corData ${outMatrix} \
                    --corMethod spearman \
                    --whatToPlot heatmap \
                    --plotFile ${outpath}/chipCorrelation/allVSall/allChipCorrelations.pdf \
                    --skipZeros \
                    --plotTitle "Correlation heatmap between ChIP samples" \
                    --plotFileFormat pdf \
                    --colorMap viridis \
                    --outFileCorMatrix ${outpath}/chipCorrelation/allVSall/allChipCorrelation_matrix.tsv 

    plotCorrelation --corData ${outMatrix} \
                    --corMethod spearman \
                    --whatToPlot scatterplot \
                    --plotFile ${outpath}/chipCorrelation/allVSall/allChipScatterPlot.pdf \
                    --skipZeros \
                    --plotTitle "Correlation scatterplot between ChIP samples" \
                    --plotFileFormat pdf
fi
