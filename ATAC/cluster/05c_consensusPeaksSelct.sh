#!/bin/bash
# -*- ENCODING: UTF-8 -*-

##===============================================================================
## SLURM VARIABLES
#SBATCH --job-name=consensusPeaks
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=0-01:00:00
#SBATCH -p short
#SBATCH -o /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.out  
#SBATCH -e /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.err 


# HOW TO RUN ME
#sbatch /home/jmendietaes/programas/pipelines/ATAC/cluster/05c_consensusPeaksSelct.sh \
#/home/jmendietaes/data/2021/ATAC/allProcessed \

# OBJECTIVE
# make consensus peaks from a selected set of peak files we have linked


# path where we have the folder structure for chip analysis 
# (inside we have ${basePath}/bamfiles/valid/)
basePath=$1
#basePath="/home/jmendietaes/data/2021/chip/allProcessed"


# Important paremeter to modify
# Set to lowercase yes to use merge of IgG from different cells allCell_IgG.sort.r..
# as control in cells with no control. Otherwise to No
# This Igg file has to be in the bams folder inside of a folder 
# called cellMergeIgG
useMergeIgG="yes"
# Set to lowercase yes to get consensus peaks of all ChIPs and cells
doAllMerge="yes"


# path for the location of the pipeline scripts
scriptsPath="/home/jmendietaes/programas/pipelines"
# species shortcut for MACS
species="mm"
speciesGenome="mm10"

# Where to look for bam files and where to store output tree
bamsPath="${basePath}/bamfiles/valid/07_LSKdays"
outpath=${basePath}"/furtherAnalysis/07_LSKdays"
inpath=${basePath}"/furtherAnalysis/02_firstATAC"

# never filter out _IgG in here
allbams=$(find ${bamsPath}/*bam -printf "${bamsPath}/%f\n" | \
            tr '\n' ' ')
allLabels=`for i in $allbams; do basename ${i} | cut -d '.' -f 1 | cut -d '_' -f 1,2,3; done | tr '\n' ' '`


## load modules
export PATH="/home/jmendietaes/programas/miniconda3/envs/Renv/bin:$PATH"
export PATH="/home/jmendietaes/programas/miniconda3/bin:$PATH"

module load SAMtools/1.12-GCC-10.2.0
module load MACS2/2.2.7.1-foss-2018b-Python-3.6.6
module load BEDTools/2.27.1-foss-2018b
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



#########################
# CONSENSUS PEAKS ANALYSIS: Consensus peaks across samples, create boolean 
#     filtering file, SAF file for featureCounts and UpSetR plot for intersection
#########################

## We first do it in the whole dataset as a trial, but in the future ill do it only between the biological replicates
if [ ! -e ${outpath}/peakCalling/MACS2/consensusPeaks ]; then
	mkdir -p ${outpath}/peakCalling/MACS2/consensusPeaks
fi

if [[ ${doAllMerge} == "yes" ]]; then

    echo -e "Starting consensus peak analysis -------------------------------------\n"


    chip="allmerged"
    for peaktype in broadPeak; do
        # select the columns to peak in each peak calling case
        if [ ${peaktype} == "narrowPeak" ]; then
            mergecols=`seq 2 10 | tr '\n' ','`
            expandparam='--is_narrow_peak'
        elif [ ${peaktype} == "broadPeak" ]; then
            mergecols=`seq 2 9 | tr '\n' ','`
            expandparam=''
        fi

        # get all non empty peak files (avoid controls)
        bamsShort=$(\
                for filename in ${allbams}; do 
                    mapLib=$(basename ${filename})
                    mapLib=(${mapLib//\./ }); 
                    mapLib=${mapLib[0]};
                    mapLib=(${mapLib//_/ }); 
                    echo ${mapLib[@]:0:2} | tr ' ' '_'; done | sort | uniq)
        
        peakFiles=$(for f in $bamsShort; do find -L ${inpath}/peakCalling/MACS2/peaks/${f}*${peaktype} \
                        -maxdepth 1  -type f ! -size 0 ; done)
        fileLabels=$(for f in $peakFiles; do echo ${f##*/} |sed "s/_peaks.${peaktype}//g"; done | tr '\n' ',')

        # check if the file exists or it was created with a previous peaks version 
        prefix="${chip}_${peaktype}_consensusPeaks"
        RoutPlot=${outpath}/peakCalling/MACS2/consensusPeaks/${prefix}.boolean.intersect.plot.pdf
        fileNotExistOrOlder "${RoutPlot}" "${peakFiles}"
        # this outputs analyse as yes or no in lowercase
        if [[ ${analyse} == "yes" ]]; then

            sort -T '.' -k1,1 -k2,2n ${peakFiles} \
                | mergeBed -c $mergecols -o collapse > ${outpath}/peakCalling/MACS2/consensusPeaks/${prefix}.txt

            python ${scriptsPath}/ChIP/cluster/02_NR_macs2_merged_expand.py ${outpath}/peakCalling/MACS2/consensusPeaks/${prefix}.txt \
                ${fileLabels} \
                ${outpath}/peakCalling/MACS2/consensusPeaks/${prefix}.boolean.txt \
                $expandparam

            consensusPeakBed=${outpath}/peakCalling/MACS2/consensusPeaks/${prefix}.bed
            awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print $1, $2, $3, $4, "0", "+" }' \
                ${outpath}/peakCalling/MACS2/consensusPeaks/${prefix}.boolean.txt > ${consensusPeakBed}
            echo -e "GeneID\tChr\tStart\tEnd\tStrand" > ${outpath}/peakCalling/MACS2/consensusPeaks/${prefix}.saf
            awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print $4, $1, $2, $3,  "+" }' \
                ${outpath}/peakCalling/MACS2/consensusPeaks/${prefix}.boolean.txt >> ${outpath}/peakCalling/MACS2/consensusPeaks/${prefix}.saf
            Rscript ${scriptsPath}/ChIP/cluster/02_NR_plot_peak_intersect.r \
                -i ${outpath}/peakCalling/MACS2/consensusPeaks/${prefix}.boolean.intersect.txt \
                -o ${RoutPlot}
        fi
        
    done

    echo -e "consensus peak analysis - Finished ------------------------------\n"
fi


############################################################
# Annotate consensus peaks with HOMER, and add annotation to boolean output file
############################################################

if [ ! -e ${outpath}/HOMER/consensusPeaks ]; then
	mkdir -p ${outpath}/HOMER/consensusPeaks
fi

if [[ ${doAllMerge} == "yes" ]]; then
    echo -e "Starting consensus peak annotations ------------------------------\n"


    chip="allmerged"
    for peaktype in broadPeak; do

        prefix="${chip}_${peaktype}_consensusPeaks"
        consensusPeakBed=${outpath}/peakCalling/MACS2/consensusPeaks/${prefix}.bed

        # check if the file exists or it was created with a previous peaks version 
        boolAnotMatr=${outpath}/HOMER/consensusPeaks/${prefix}.boolean.annotatePeaks.txt
        fileNotExistOrOlder "${boolAnotMatr}" "${consensusPeakBed}"
        # this outputs analyse as yes or no in lowercase
        if [[ ${analyse} == "yes" ]]; then

            annotatePeaks.pl \
                    ${consensusPeakBed} \
                    ${speciesGenome} \
                    -gid \
                    -cpu ${SLURM_CPUS_PER_TASK} \
                    -annStats ${outpath}/HOMER/consensusPeaks/${prefix}.annotateStats.txt \
                    > ${outpath}/HOMER/consensusPeaks/${prefix}.annotatePeaks.txt

            cut -f2- ${outpath}/HOMER/consensusPeaks/${prefix}.annotatePeaks.txt | \
                awk 'NR==1; NR > 1 {print $0 | "sort -T '.' -k1,1 -k2,2n"}' | \
                cut -f6- > ${outpath}/HOMER/consensusPeaks/tmp.txt
            paste ${outpath}/peakCalling/MACS2/consensusPeaks/${prefix}.boolean.txt \
                ${outpath}/HOMER/consensusPeaks/tmp.txt > ${boolAnotMatr}
        fi
    done

    echo -e "consensus peak annotations - Finished ---------------------\n"
fi


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
    for peaktype in broadPeak; do
        prefix="${chip}_${peaktype}_consensusPeaks"
        # get all non empty peak files (avoid controls)
        bamsShort=$(\
                for filename in ${allbams}; do 
                    mapLib=$(basename ${filename})
                    mapLib=(${mapLib//\./ }); 
                    mapLib=${mapLib[0]};
                    mapLib=(${mapLib//_/ }); 
                    echo ${mapLib[@]:0:2} | tr ' ' '_'; done | sort | uniq)
        
        peakFiles=$(for f in $bamsShort; do find -L ${inpath}/peakCalling/MACS2/peaks/${f}*${peaktype} \
                        -maxdepth 1  -type f ! -size 0 | { grep -v -e "_input" -v -e "_IgG" || :; }; done)
        fileLabels=$(for f in $peakFiles; do echo ${f##*/} |sed "s/_peaks.${peaktype}//g"; done)
        # get a list with the bam files in the same order as peak files
        bamfiles=$(for f in ${fileLabels}; do find -L ${bamsPath}/${f}*bam -printf "${bamsPath}/%f "; 
                        done | tr '\n' ' ')


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


###########################################################
# Read count to CPM
###########################################################

if [[ ${doAllMerge} == "yes" ]]; then
    echo -e "Starting consensus CPM -----------------------\n"
    chip="allmerged"
    for peaktype in broadPeak; do
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



exit 0
