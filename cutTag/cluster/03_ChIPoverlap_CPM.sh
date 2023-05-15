#!/bin/bash
# -*- ENCODING: UTF-8 -*-

##===============================================================================
## SLURM VARIABLES
#SBATCH --job-name=cutTag_CPM
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
#SBATCH --time=10:00:00
#SBATCH -p short
#SBATCH -o /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.out  
#SBATCH -e /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.err 


# HOW TO RUN ME
#sbatch /home/jmendietaes/programas/pipelines/cutTag/cluster/03_ChIPoverlap_CPM.sh \
#/home/jmendietaes/data/2021/cutTag/allProcessed \
# OBJECTIVE
# get cutTag CPM at location overlaping ChIP peaks for same factor
# HERE WE DONT USE SPIKE-IN INFO (MAINLY BECAUSE IS NOT WORKING FOR US)

## How to adapt epic peaks
# Prepare Epic files
# filename=DM_K9me3_epic.out
# filenameOut=DM_K9me3_epic_peaks.broadPeak

# nPeak=$(wc -l ${filename}  | awk '{print $1}') ; 
# #echo -e "peakName" > temp1.txt ; 
# for i in `seq 1 $(($nPeak-1))`; do echo Epic_peak_${i}; done > temp1.txt

# grep -v "^#" ${filename} | awk '{print $5"\t.\t"$10"\t"$4"\t"$9}' > temp2.txt

# grep -v "^#" ${filename} | awk '{print $1"\t"$2"\t"$3}' | \
#     paste - temp1.txt temp2.txt > ${filenameOut}

# rm temp1.txt; rm temp2.txt


# path where we have the folder structure for cutTag analysis 
# (inside we have ${basePath}/bamfiles/valid/)
basePath=$1
#basePath="/home/jmendietaes/data/2021/cutTag/allProcessed"

# Path to ChIP peaks
peakPath="/home/jmendietaes/data/2021/chip/allProcessed/furtherAnalysis/08_projectRestart/peakCalling/MACS2/peaks"
# path for the location of the pipeline scripts
scriptsPath="/home/jmendietaes/programas/PhD"

# Where to look for bam files and where to store output tree
bamsPath="${basePath}/bamfiles/valid/03_properAnalysis"
outpath=${basePath}"/furtherAnalysis/03_properAnalysis"


# Get list of focus bam files
allbams=$(find ${bamsPath}/*bam -printf "${bamsPath}/%f\n" | \
            tr '\n' ' ')


allLabels=`for i in $allbams; do basename ${i} | cut -d '.' -f 1 | \
            cut -d '_' -f 1,2,3; done | tr '\n' ' '`
#SLURM_CPUS_PER_TASK=6

# First get a list of all ChIP
chips=$(for a in ${allLabels}; do 
    epiName=(${a//_/ }); epiName=${epiName[1]}; 
    epiName=(${epiName//-/ }); epiName=${epiName[0]};
    echo ${epiName}; done | sort | uniq)


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


#########################
# CONSENSUS Folder: Define or recycle the output path for these specific
# consensus
#########################

#### Get the output path for these consensus
if [ ! -e ${outpath}/specificConsensus ]; then
    mkdir -p ${outpath}/specificConsensus
fi

for chip in ${chips}; do
    echo "############  ${chip}  ############"

    peakLabels=$(echo "${allLabels}" | tr ' ' '\n' | grep $chip)
    
    
    peakCells=$(for a in ${peakLabels}; do 
            cell=(${a//_/ }); cell=${cell[0]}; 
            cell=(${cell//-/ }); cell=${cell[0]};
            echo ${cell}; done | sort | uniq)


    for cell in ${peakCells}; do 
        echo $cell

        peakLabelsCell=$(echo $peakLabels | tr ' ' '\n' | \
                    { grep "${cell}-"* || :; } )
        echo ${peakLabelsCell}
        # get a list with the bam files in the same order as peak labels
        bamfiles=$(for f in ${peakLabelsCell}; do find -L ${bamsPath}/${f}*bam \
                    -printf "${bamsPath}/%f "; 
                        done | tr '\n' ' ')

        peakLabelsTxt=$(echo $peakLabelsCell | tr ' ' ':' )

        # if empty list, pass
        nchar=$(echo $peakLabelsCell | wc -w)
        if [[ "$nchar" -ne 0 ]] ; then

            unset -v latest
            for file in "${outpath}/specificConsensus/"consensus-${cell}-${chip}_*; do
            [[ $file -nt $latest ]] && latest=$file
            done
            echo $latest

            # if empty folder define first output name
            nchar=$(echo $latest | wc -w)
            if [[ "$nchar" -eq 0 ]] ; then
                consensusOut="${outpath}/specificConsensus/consensus-${cell}-${chip}_1";
                mkdir ${consensusOut}
                echo ${peakLabelsTxt} > ${consensusOut}/combinedPeakFiles.txt
            else
                restarted="no"
                # If there were previous files maybe we had the same comparison
                #once (we can try to reuse what is newer than the bam or peak files)
                for file in "${outpath}/specificConsensus/"consensus-${cell}-${chip}_*; do
                    checkLbl=`sed "1q;d" ${file}/combinedPeakFiles.txt`
                    basefi=$(basename ${file})
                    if [[ ${peakLabelsTxt} == ${checkLbl} ]]; then 
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
                    consensusOut="${outpath}/specificConsensus/consensus-${cell}-${chip}_${lastN}"
                    mkdir ${consensusOut}
                    echo ${peakLabelsTxt} > ${consensusOut}/combinedPeakFiles.txt
                fi
            fi


            echo "Output consensus folder: "
            echo ${consensusOut}




            #########################
            # CONSENSUS PEAKS ANALYSIS: Get ChIP peak coordiantes and store them as 
            #   saf
            #########################

            ## We first do it in the whole dataset as a trial, but in the future ill do it only between the biological replicates
            if [ ! -e ${consensusOut}/consensusPeaks ]; then
                mkdir -p ${consensusOut}/consensusPeaks
            fi


            echo -e "Starting ChIP peak retrieval -------------------------------------\n"


            for peaktype in broadPeak; do
                peakP=$(find ${peakPath}/${cell}_*${chip}*${peaktype} \
                    -printf "${peakPath}/%f\n" | \
                    { grep -e "${chip}_" -e "${chip}-" || :; } | tr '\n' ' ')
                peakName=$(basename ${peakP}).saf

                # check if the file exists or it was created with a previous peaks version 
                fileNotExistOrOlder "${consensusOut}/consensusPeaks/${peakName}" \
                                    "${peakP}"
                # this outputs analyse as yes or no in lowercase
                if [[ ${analyse} == "yes" ]]; then
                    echo -e "Geneid\tChr\tStart\tEnd\tStrand" > \
                        ${consensusOut}/consensusPeaks/${peakName}
                    awk '{print $4"\t"$1"\t"$2"\t"$3"\t+"}' ${peakP} >> \
                        ${consensusOut}/consensusPeaks/${peakName}
                fi

            done


            echo -e "ChIP peak retrieval - Finished ------------------------------\n"



            ###########################################################
            # Count reads in consensus peaks with featureCounts
            ##########################################################

            featureCpath="${consensusOut}/consensusPeaks/featureCounts"
            if [ ! -e ${featureCpath} ]; then
                mkdir -p ${featureCpath}
            fi
            cd ${featureCpath}

            echo -e "Starting consensus featureCounts -----------------------\n"

            for peaktype in broadPeak; do
                prefix="${cell}_${chip}_${peaktype}_consensusPeaks"
                
                # check if the file exists or it was created with a previous peaks version 
                featureOut=${featureCpath}/${prefix}.featureCounts.txt
                fileNotExistOrOlder "${featureOut}" "${peakFiles}"
                # this outputs analyse as yes or no in lowercase
                if [[ ${analyse} == "yes" ]]; then
                    # this only for the consensus between replicates, here no sense
                    head ${consensusOut}/consensusPeaks/${peakName};
                    featureCounts \
                            -F SAF \
                            -O \
                            --fracOverlap 0.2 \
                            -T ${SLURM_CPUS_PER_TASK} \
                            -p --donotsort \
                            -a ${consensusOut}/consensusPeaks/${peakName} \
                            -o ${featureOut} \
                            ${bamfiles}
                fi

            done
            echo -e "Consensus featureCounts - Finished ----------------------\n"



            ###########################################################
            # Read count to CPM
            ###########################################################

            echo -e "Starting consensus CPM -----------------------\n"
            for peaktype in broadPeak; do
                prefix="${cell}_${chip}_${peaktype}_consensusPeaks"
                
                # featureCount files start with this columns before the bam read counts
                #Geneid	Chr	Start	End	Strand	Length
                prevFields=6
                lengthCol=5

                # set the output file path and copy the content of original
                featureCPM=${featureCpath}/${prefix}.featureCounts.CPM.txt

                echo -e "FileName\tSampleName\tCellType\tStatus" > ${featureCpath}/sampleInfo.txt
                for bam in ${bamfiles}; do
                    sname=`basename $bam  | sed 's/.sort.rmUnM.q30.rmchr.Tn5.bam//g'`
                    #sname=$bam
                    # cell is always the same here
                    #cell=(${sname//-/ }) ; cell=${cell[0]}
                    status=(${sname//_/ }) ; status=${status[0]}; 
                    status=(${status//-/ }); 
                    status=$(echo ${status[@]:1} | tr ' ' '-')

                    echo -e ${bam}"\t"${sname}"\t"${cell}"\t"${status} >> ${featureCpath}/sampleInfo.txt
                done

                # check if the file exists or it was created with a previous featureCounts version 
                fileNotExistOrOlder "${featureCPM}" "${featureCpath}/${prefix}.featureCounts.txt"
                # this outputs analyse as yes or no in lowercase
                if [[ ${analyse} == "yes" ]]; then
                    # this only for the consensus between replicates, here no sense
                    Rscript ${scriptsPath}/cutTag/cluster/02_NR_metric_CPMfromFeatureCounts.r \
                                -i ${featureCpath}/${prefix}.featureCounts.txt \
                                -o ${featureCPM} \
                                -s ${featureCpath}/sampleInfo.txt \
                                -l ${lengthCol} \
                                -d ${prevFields}
                                    

                fi

            done
            echo -e "Consensus CPM - Finished ----------------------\n"


            ############################################
            # Adding back ChIP peak calling info
            ############################################
            echo -e "Starting ChIP reinfo ---------------------------\n"

            # for now we only work with narrow peaks
            for peaktype in broadPeak; do


                prefix="${cell}_${chip}_${peaktype}_consensusPeaks"
                outfile=${featureCpath}/${cell}_${chip}_${peaktype}_cutTag.overlap.txt

                echo -e "score\tsignalValue\tpValue\tqValue" > ${featureCpath}/tmp.txt
                awk '{print $5"\t"$7"\t"$8"\t"$9}' ${peakP} >> ${featureCpath}/tmp.txt

                awk {'first = $1; $1=""; print $0'} ${featureCPM} | \
                    sed 's/^ //g' | sed 's/ /\t/g'> ${featureCpath}/tmp2.txt
                headerNames=$(head -n 1 ${featureCpath}/tmp2.txt)
                for h in ${headerNames}; do
                    sed -i "s/${h}/${h}.cpm/g" ${featureCpath}/tmp2.txt
                done
                # Both files are in the same order
                tail -n +2 ${featureCpath}/${prefix}.featureCounts.txt | \
                    sed "s/sort.rmUnM.q30.rmchr.Tn5.bam/rcount/g" | \
                    sed "s/$(echo $bamsPath/ | sed 's_/_\\/_g')//g" | \
                    paste - \
                    ${featureCpath}/tmp2.txt ${featureCpath}/tmp.txt > ${outfile}


            done
        fi
    done
done
echo -e "END --------------------------------------------------"

