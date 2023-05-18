#!/bin/bash
# -*- ENCODING: UTF-8 -*-

##===============================================================================
## SLURM VARIABLES
#SBATCH --job-name=countIgG
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH -p short
#SBATCH -o /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.out  
#SBATCH -e /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.err 


# HOW TO RUN ME
#sbatch /home/jmendietaes/programas/pipelines/ChIP/cluster/05b_countReadsIgG.sh \
#/home/jmendietaes/data/2021/chip/allProcessed \

# OBJECTIVE
# Get read counts and CPM of IgG in consensus peak coordinates
# overwrites previous results to add the IgG data

# Important paremeter to modify
# Set to lowercase yes to use merge of IgG from different cells allCell_IgG.sort.r..
# as control in cells with no control. Otherwise to No
# This Igg file has to be in the bams folder inside of a folder 
# called cellMergeIgG
useMergeIgG="yes"
# Set to lowercase yes to get consensus peaks of all ChIPs and cells
doAllMerge="yes"

# path where we have the folder structure for chip analysis 
# (inside we have ${basePath}/bamfiles/valid/)
basePath=$1
#basePath="/home/jmendietaes/data/2021/chip/allProcessed"
# Location of the analysis header files downloaded from 
# https://github.com/nf-core/chipseq/tree/master/assets/multiqc
#extraFilePath=$2
#extraFilePath="/home/jmendietaes/data/2021/chip/analysisFiles"


# path for the location of the pipeline scripts
scriptsPath="/home/jmendietaes/programas/pipelines"
# extend variables
bamsPath="${basePath}/bamfiles/valid/subsampled_noIgG"
outpath=${basePath}"/furtherAnalysis/subsampled_noIgG"

# never filter out _IgG in here
allbams=$(find ${bamsPath}/*bam -printf "${bamsPath}/%f\n" | \
            tr '\n' ' ')
allLabels=`for i in $allbams; do basename ${i} | cut -d '.' -f 1 | cut -d '_' -f 1,2,3; done | tr '\n' ' '`
#SLURM_CPUS_PER_TASK=6


## load modules
export PATH="/home/jmendietaes/programas/miniconda3/envs/Renv/bin:$PATH"
export PATH="/home/jmendietaes/programas/miniconda3/bin:$PATH"

module load SAMtools/1.12-GCC-10.2.0
module load MACS2/2.2.7.1-foss-2018b-Python-3.6.6
module load BEDTools/2.27.1-foss-2018b
module load Homer/4.10-foss-2018b
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

# function to add their respective IgG or control files to a 
# list of bamfiles we have (to get the reads by featureCounts)
# first input is bamfiles variable with the bams we are looking at
# second input is allbams variable with all bam files we have
# third input is useMergeIgG
# 4th input is bamsPath with the path to all bamfiles
addControls () {
# now we want to also include IgG and input files coverage
    # so we look to all the cells we have
    checkCells=()
    for bam in $1; do
        cell=$(basename $bam | cut -d '.' -f 1 | cut -d '_' -f 1)
        if [[ ! " ${checkCells[*]} " =~ " ${cell} " ]]; then
            checkCells+=("${cell}")
        fi
    done
    # and add the corresponding controls
    cellControls=()
    for cell in ${checkCells[*]}; do
        cellC=`echo ${2} | tr ' ' '\n' | grep "${cell}_" | \
                        { grep -e "_input" -e "_IgG" || :; }`

        # IMPORTANT
        # if we have no control for this cell and stated useMergeIgG to Yes we use the merge
        if [[ $3 == "yes" ]] && [[ $cellC == "" ]]; then
            cellC=`find ${4}/cellMergeIgG/allCell_*bam`
        fi
        if [[ ! " ${cellControls[*]} " =~ " ${cellC} " ]]; then
            cellControls+=("${cellC}")
        fi
    done

    # rename bamfiles
    bamfiles=$(printf "%s "  $1 "${cellControls[@]}")
    bamfiles=$(for ba in ${bamfiles}; do echo $ba; done | tr '\n' ' ')
}


# first make sure input and output folder do not contain names
# that can result in an issue
if [[ $bamsPath == *"_IgG"* ]]; then
    echo "bamsPath contains key substring _IgG";
    exit 1;
fi
if [[ $outpath == *"_IgG"* ]]; then
    echo "outpath contains key substring _IgG";
    exit 1;
fi

if [[ $bamsPath == *"_input"* ]]; then
    echo "bamsPath contains key substring _input"
    exit 1;
fi
if [[ $outpath == *"_input"* ]]; then
    echo "outpath contains key substring _input"
    exit 1;
fi


################################# START #######################
###########################################################
# Count reads in consensus peaks with featureCounts
##########################################################

featureCpath=${outpath}/peakCalling/MACS2/consensusPeaks/featureCounts
if [ ! -e ${featureCpath} ]; then
	mkdir -p ${featureCpath}
fi
cd ${featureCpath}

if [[ ${doAllMerge} == "yes" ]]; then
    echo -e "Starting consensus featureCounts -----------------------\n"

    chip="allmerged"
    for peaktype in narrowPeak broadPeak; do
        prefix="${chip}_${peaktype}_consensusPeaks"
        peakFiles=$(find -L ${outpath}/peakCalling/MACS2/peaks/*.${peaktype} -maxdepth 1  -type f ! -size 0 | \
                    { grep -v -e "_input" -v -e "_IgG" || :; } )
        fileLabels=$(for f in $peakFiles; do echo ${f##*/} |sed "s/_peaks.${peaktype}//g"; done)
        # get a list with the bam files in the same order as peak files
        bamfiles=$(for f in ${fileLabels}; do find -L ${bamsPath}/${f}*bam -printf "${bamsPath}/%f "; 
                        done | tr '\n' ' ')
        # now we want to also include IgG files in bamfiles variable
        addControls "${bamfiles}" "${allbams}" "${useMergeIgG}" "${bamsPath}"


        # check if the file exists or it was created with a previous peaks version 
        featureOut=${featureCpath}/${prefix}.featureCounts.txt
        fileNotExistOrOlder "${featureOut}" "${peakFiles}"
        # this outputs analyse as yes or no in lowercase
        if [[ ${analyse} == "yes" ]]; then
            # this only for the consensus between replicates, here no sense
            featureCounts \
                    -F SAF \
                    -O \
                    --fracOverlap 0.2 \
                    -T ${SLURM_CPUS_PER_TASK} \
                    -p --donotsort \
                    -a ${outpath}/peakCalling/MACS2/consensusPeaks/${prefix}.saf \
                    -o ${featureOut} \
                    ${bamfiles}
        fi
    done
    echo -e "Consensus featureCounts - Finished ----------------------\n"
fi

## Same by comparing samples of same chip

featureCPairpath=${outpath}/peakCalling/MACS2/consensusPeaks/bySameChip/featureCounts
if [ ! -e ${featureCPairpath} ]; then
	mkdir -p ${featureCPairpath}
fi
cd ${featureCPairpath}

echo -e "Starting consensus same-chip featureCounts -----------------------\n"

# get all analysed ChIP
sameChipCons=${outpath}/peakCalling/MACS2/consensusPeaks/bySameChip
more1Chip=$(\
    for filename in ${sameChipCons}/*bed; do 
        filename=$(basename ${filename})
        mapLib=(${filename//_/ }); 
        echo ${mapLib[0]}; done | \
        sort | uniq)

for chip in ${more1Chip}; do
    for peaktype in narrowPeak broadPeak; do
        prefix="${chip}_${peaktype}_consensusPeaks"
        peakFiles=$(find -L ${outpath}/peakCalling/MACS2/peaks/*.${peaktype} -maxdepth 1  -type f ! -size 0 | \
                    { grep "${chip}" || :; } )
        fileLabels=$(for f in $peakFiles; do echo ${f##*/} |sed "s/_peaks.${peaktype}//g"; done)
        # get a list with the bam files in the same order as peak files
        bamfiles=$(for f in ${fileLabels}; do find -L ${bamsPath}/${f}*bam -printf "${bamsPath}/%f "; 
                        done | tr '\n' ' ')

        # now we want to also include IgG files in bamfiles variable
        addControls "${bamfiles}" "${allbams}" "${useMergeIgG}" "${bamsPath}"

        # check if the file exists or it was created with a previous peaks version 
        featureOut=${featureCPairpath}/${prefix}.featureCounts.txt
        fileNotExistOrOlder "${featureOut}" "${peakFiles}"
        # this outputs analyse as yes or no in lowercase
        if [[ ${analyse} == "yes" ]]; then
            # this only for the consensus between replicates, here no sense
            featureCounts \
                    -F SAF \
                    -O \
                    --fracOverlap 0.2 \
                    -T ${SLURM_CPUS_PER_TASK} \
                    -p --donotsort \
                    -a ${outpath}/peakCalling/MACS2/consensusPeaks/bySameChip/${prefix}.saf \
                    -o ${featureOut} \
                    ${bamfiles} 
        fi
    done
done

echo -e "Consensus same-chip featureCounts - Finished ----------------------\n"

###########################################################
# Read count to CPM
###########################################################

if [[ ${doAllMerge} == "yes" ]]; then
    echo -e "Starting consensus CPM -----------------------\n"
    chip="allmerged"
    for peaktype in narrowPeak broadPeak; do
        prefix="${chip}_${peaktype}_consensusPeaks"
        
        bamfiles=$(head -n 2 ${featureCpath}/${prefix}.featureCounts.txt | tail -n 1 |\
                    tr '\t' '\n' | grep bam$)
        nbamfiles=$(echo $bamfiles | wc -w)
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
                        -d ${prevFields} \
                         

        fi
    done
    echo -e "Consensus CPM - Finished ----------------------\n"
fi

## Same by comparing samples of same chip
echo -e "Starting same-chip CPM -----------------------\n"

for chip in ${more1Chip}; do
    for peaktype in narrowPeak broadPeak; do
        prefix="${chip}_${peaktype}_consensusPeaks"
        
        bamfiles=$(head -n 2 ${featureCPairpath}/${prefix}.featureCounts.txt | tail -n 1 |\
                    tr '\t' '\n' | grep bam$)
        nbamfiles=$(echo $bamfiles | wc -w)
        # featureCount files start with this columns before the bam read counts
        #Geneid	Chr	Start	End	Strand	Length
        prevFields=6
        lengthCol=5

        # set the output file path and copy the content of original
        featureCPM=${featureCPairpath}/${prefix}.featureCounts.CPM.txt

        echo -e "FileName\tSampleName\tCellType\tStatus" > ${featureCPairpath}/sampleInfo.txt
        for bam in ${bamfiles}; do
            sname=$(basename $bam | cut -d '.' -f 1 | cut -d '_' -f 1,2,3)
            #sname=$bam
            cell=(${sname//_/ }) ; cell=${cell[0]}
            status=(${sname//_/ }) ; status=${status[1]}; status=(${status//-/ }); status=${status[0]}
            echo -e ${bam}"\t"${sname}"\t"${cell}"\t"${status} >> ${featureCPairpath}/sampleInfo.txt
        done

        # check if the file exists or it was created with a previous featureCounts version 
        fileNotExistOrOlder "${featureCPM}" "${featureCPairpath}/${prefix}.featureCounts.txt"
        # this outputs analyse as yes or no in lowercase
        if [[ ${analyse} == "yes" ]]; then
            # this only for the consensus between replicates, here no sense
            Rscript ${scriptsPath}/ChIP/cluster/02_NR_metric_CPMfromFeatureCounts.r \
                    -i ${featureCPairpath}/${prefix}.featureCounts.txt \
                    -o ${featureCPM} \
                    -s ${featureCPairpath}/sampleInfo.txt \
                    -l ${lengthCol} \
                    -d ${prevFields}

        fi
    done
done
echo -e "Same-chip CPM - Finished ----------------------\n"


###########################################################
# Make final table with all data
##########################################################

if [[ ${doAllMerge} == "yes" ]]; then
    echo -e "Starting final merge table -----------------------\n"

    # takes almost no time, so will repeat it always
    python ${scriptsPath}/ChIP/cluster/02_NR_gatherAllInTable.py \
                ${outpath}/HOMER/consensusPeaks \
                ${outpath}/peakCalling/MACS2/consensusPeaks/featureCounts \
                ${outpath}/peakCalling/MACS2/consensusPeaks/annotWholeTable

    echo -e "Final merge table- Finished ----------------------\n"
fi

## Same by comparing samples of same chip
echo -e "Starting same-chip final merge table -----------------------\n"

# takes almost no time, so will repeat it always
python ${scriptsPath}/ChIP/cluster/02_NR_gatherAllInTable.py \
            ${outpath}/HOMER/consensusPeaks/bySameChip \
            ${outpath}/peakCalling/MACS2/consensusPeaks/bySameChip/featureCounts \
            ${outpath}/peakCalling/MACS2/consensusPeaks/bySameChip/annotWholeTable

echo -e "Same-chip final merge table - Finished ----------------------\n"

