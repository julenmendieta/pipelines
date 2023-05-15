#!/bin/bash
# -*- ENCODING: UTF-8 -*-

##===============================================================================
## SLURM VARIABLES
#SBATCH --job-name=specificConsensus
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
#SBATCH --time=10:00:00
#SBATCH -p short
#SBATCH -o /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.out  
#SBATCH -e /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.err 


# HOW TO RUN ME
#sbatch /home/jmendietaes/programas/pipelines/ChIP/cluster/07a_specificConsensus_bigwigTMM.sh \
#/home/jmendietaes/data/2021/chip/allProcessed \

# OBJECTIVE
# get scaling factors for group of cells and ChIPs we want to compare in the
#future. Then generate scaled bigwig files


# path where we have the folder structure for chip analysis 
# (inside we have ${basePath}/bamfiles/valid/)
basePath=$1
#basePath="/home/jmendietaes/data/2021/chip/allProcessed"
# Location of the analysis header files downloaded from 
# https://github.com/nf-core/chipseq/tree/master/assets/multiqc
#extraFilePath=$2
#extraFilePath="/home/jmendietaes/data/2021/chip/analysisFiles"


# path for the location of the pipeline scripts
scriptsPath="/home/jmendietaes/programas/PhD"

# Where to look for bam files and where to store output tree
bamsPath="${basePath}/bamfiles/valid/subsampled_noIgG"
outpath=${basePath}"/furtherAnalysis/subsampled_noIgG"

# Space separated list of files for the consensus and future analysis
# With unique begining of ID is fine (2-3 first sections in between _)
# E.j.: "MEP_Brd9-merged-sub105350855 GMP_Smarcb1_ChIP11"
# Actually if fails if 4 sections are used. Ej Bcell_Brd9_ChIP9_S11
# instead of Bcell_Brd9_ChIP9
focusFiles="Bcell_Brd9_ChIP9 \
Bcell_Kmt2a_ChIP12 \
Bcell_Kmt2d_ChIP12 \
Ery_Brd9-merged \
Ery_Kmt2a_ChIP11 \
Ery_Kmt2d-merged \
Ery_Smarcb1_ChIP9 \
GMP_Brd9-merged \
GMP_Kmt2a_ChIP11 \
GMP_Kmt2d_ChIP11 \
GMP_Smarcb1_ChIP11 \
MEP_Brd9-merged \
MEP_Kmt2a_ChIP11 \
MEP_Kmt2d_ChIP11 \
MEP_Smarcb1_ChIP12 \
Mono_Brd9-merged \
Mono_Kmt2a-merged \
Mono_Kmt2d-merged \
Mono_Smarcb1-merged"



# Get list of focus bam files
allbams=$(find ${bamsPath}/*bam -printf "${bamsPath}/%f\n" | \
            tr '\n' ' ')
focusGrep=$(echo $focusFiles | sed 's/ /\\|/g')
allbams=$(echo $allbams | tr ' ' '\n' | grep $focusGrep)

allLabels=`for i in $allbams; do basename ${i} | cut -d '.' -f 1 | cut -d '_' -f 1,2,3; done | tr '\n' ' '`
allLabelsTxt=$(echo $allLabels | tr ' ' ':')
#SLURM_CPUS_PER_TASK=6


## load modules
export PATH="/home/jmendietaes/programas/miniconda3/envs/Renv/bin:$PATH"
export PATH="/home/jmendietaes/programas/miniconda3/bin:$PATH"

module load SAMtools/1.12-GCC-10.2.0
module load BEDTools/2.27.1-foss-2018b

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

# first make sure input and output folder do not contain names
# that can result in an issue
if [[ $bamsPath == *"_IgG"* ]]; then
    echo "bamsPath contains key substring: _IgG";
    exit 1;
fi
if [[ $outpath == *"_IgG"* ]]; then
    echo "outpath contains key substring: _IgG";
    exit 1;
fi

if [[ $bamsPath == *"_input"* ]]; then
    echo "bamsPath contains key substring: _input"
    exit 1;
fi
if [[ $outpath == *"_input"* ]]; then
    echo "outpath contains key substring: _input"
    exit 1;
fi

#########################
# CONSENSUS Folder: Define or recycle the output path for these specific
# consensus
#########################

#### Get the output path for these consensus
if [ ! -e ${outpath}/specificConsensus ]; then
	mkdir -p ${outpath}/specificConsensus
fi

unset -v latest
for file in "${outpath}/specificConsensus/"consensusSet_*; do
  [[ $file -nt $latest ]] && latest=$file
done
echo $latest

# if empty folder define first output name

nchar=$(echo $latest | wc -w)
if [[ "$nchar" -eq 0 ]] ; then
    consensusOut="${outpath}/specificConsensus/consensusSet_1";
    mkdir ${consensusOut}
    echo ${allLabelsTxt} > ${consensusOut}/combinedPeakFiles.txt
else
    restarted="no"
    # If there were previous files maybe we had the same comparison
    #once (we can try to reuse what is newer than the bam or peak files)
    for file in "${outpath}/specificConsensus/"consensusSet_*; do
        checkLbl=`sed "1q;d" ${file}/combinedPeakFiles.txt`
        basefi=$(basename ${file})
        if [[ ${allLabelsTxt} == ${checkLbl} ]]; then 
            restarted="yes"
            consensusOut="${outpath}/specificConsensus/${basefi}";
            echo "Using previous consensus folder"
        fi
    done

    # If there are more consensus folders but none compare same chips
    if [ ${restarted} == "no" ]; then
        lastN=$(basename ${latest})
        lastN=(${lastN//_/ }); lastN=${lastN[1]}
        lastN=$(echo "print($lastN + 1)" | python3)
        consensusOut="${outpath}/specificConsensus/consensusSet_${lastN}"
        mkdir ${consensusOut}
        echo ${allLabelsTxt} > ${consensusOut}/combinedPeakFiles.txt
    fi
fi


echo "Output consensus folder: "
echo ${consensusOut}




#########################
# CONSENSUS PEAKS ANALYSIS: Consensus peaks across samples, create boolean 
#     filtering file, SAF file for featureCounts and UpSetR plot for intersection
#########################

## We first do it in the whole dataset as a trial, but in the future ill do it only between the biological replicates
if [ ! -e ${consensusOut}/consensusPeaks ]; then
    mkdir -p ${consensusOut}/consensusPeaks
fi


echo -e "Starting consensus peak analysis -------------------------------------\n"


chip="merged"
for peaktype in narrowPeak broadPeak; do
    # select the columns to peak in each peak calling case
    if [ ${peaktype} == "narrowPeak" ]; then
        mergecols=`seq 2 10 | tr '\n' ','`
        expandparam='--is_narrow_peak'
    elif [ ${peaktype} == "broadPeak" ]; then
        mergecols=`seq 2 9 | tr '\n' ','`
        expandparam=''
    fi

    # get all non empty peak files (avoid controls)
    peakFiles=$(find -L ${outpath}/peakCalling/MACS2/peaks/*.${peaktype} -maxdepth 1  -type f ! -size 0 | \
                { grep -e ${focusGrep} || :; } )
    fileLabels=$(for f in $peakFiles; do echo ${f##*/} |sed "s/_peaks.${peaktype}//g"; done | tr '\n' ',')

    # check if the file exists or it was created with a previous peaks version 
    prefix="${chip}_${peaktype}_consensusPeaks"
    consensusSaf=${consensusOut}/consensusPeaks/${prefix}.saf
    fileNotExistOrOlder "${consensusSaf}" "${peakFiles}"
    # this outputs analyse as yes or no in lowercase
    if [[ ${analyse} == "yes" ]]; then

        echo ${consensusSaf}

        sort -T '.' -k1,1 -k2,2n ${peakFiles} \
            | mergeBed -c $mergecols -o collapse > ${consensusOut}/consensusPeaks/${prefix}.txt

        python ${scriptsPath}/ChIP/cluster/02_NR_macs2_merged_expand.py ${consensusOut}/consensusPeaks/${prefix}.txt \
            ${fileLabels} \
            ${consensusOut}/consensusPeaks/${prefix}.boolean.txt \
            $expandparam

        consensusPeakBed=${consensusOut}/consensusPeaks/${prefix}.bed
        awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print $1, $2, $3, $4, "0", "+" }' \
            ${consensusOut}/consensusPeaks/${prefix}.boolean.txt > ${consensusPeakBed}
        echo -e "GeneID\tChr\tStart\tEnd\tStrand" > ${consensusSaf}
        awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print $4, $1, $2, $3,  "+" }' \
            ${consensusOut}/consensusPeaks/${prefix}.boolean.txt >> ${consensusSaf}

    fi
    
done

echo -e "consensus peak analysis - Finished ------------------------------\n"



###########################################################
# Count reads in consensus peaks with featureCounts
##########################################################

featureCpath="${consensusOut}/consensusPeaks/featureCounts"
if [ ! -e ${featureCpath} ]; then
    mkdir -p ${featureCpath}
fi
cd ${featureCpath}

echo -e "Starting consensus featureCounts -----------------------\n"

chip="merged"
for peaktype in narrowPeak broadPeak; do
    prefix="${chip}_${peaktype}_consensusPeaks"
    peakFiles=$(find -L ${outpath}/peakCalling/MACS2/peaks/*.${peaktype} -maxdepth 1  -type f ! -size 0 | \
                { grep -e ${focusGrep} || :; } )
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
                -a ${consensusOut}/consensusPeaks/${prefix}.saf \
                -o ${featureOut} \
                ${bamfiles}
    fi
done
echo -e "Consensus featureCounts - Finished ----------------------\n"



###########################################################
# Read count to CPM
###########################################################

echo -e "Starting consensus CPM -----------------------\n"
chip="merged"
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


############################################
# BigWig normalisation
############################################
echo -e "Starting BigWigs scaling ---------------------------\n"

# Create output dir
if [ ! -e ${consensusOut}/bigWig_TMM ]; then
    mkdir -p ${consensusOut}/bigWig_TMM
fi

# for now we only work with narrow peaks
chip="merged"
for peaktype in narrowPeak; do
    # Create output dir
    if [ ! -e ${consensusOut}/bigWig_TMM/${peaktype} ]; then
        mkdir -p ${consensusOut}/bigWig_TMM/${peaktype}
    fi

    prefix="${chip}_${peaktype}_consensusPeaks"
    scalingVals=${featureCpath}/${prefix}.featureCounts.CPM.TMMscale.txt

    # check content of eleventh line of step control file
    for fi in ${allbams}; do
        echo ${fi}
        fi=$(basename $fi)
        fi=(${fi//\.bam/ }); fi=${fi[0]}

        # get string of interest from name
        id1=(${fi//\./ }); id1=${id1[0]}
        id1=(${id1//-/ }); id1=${id1[0]}
        id1=(${id1//_/ }); id1="${id1[0]}_${id1[1]}"
        
        scalingFactor=$(grep ${id1} ${scalingVals} | cut -f 2)

        bamPath="${bamsPath}/${fi}.bam"
        bigWigOut2="${consensusOut}/bigWig_TMM/${peaktype}/${fi}.TMM.bw"

        # check if the file exists or it was created with a previous bam version 
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

exit 0