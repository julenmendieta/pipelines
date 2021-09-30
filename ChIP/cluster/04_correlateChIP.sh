#!/bin/bash
# -*- ENCODING: UTF-8 -*-

##===============================================================================
## SLURM VARIABLES
#SBATCH --job-name=bamChipCor
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH -p medium
#SBATCH -o /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.out  
#SBATCH -e /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.err 


# HOW TO RUN ME
# create file with ChIP names in bams folder
# for filename in *bam; do mapLib=(${filename//_/ }); mapLib=${mapLib[1]}; mapLib=(${mapLib//-/ }); echo ${mapLib[0]}; done | sort | uniq > chipMarcs.txt
# N=`cat chipMarcs.txt | wc -l`
#sbatch --array=1-${N} /home/jmendietaes/programas/PhD/ChIP/cluster/04_correlateChIP.sh \
#/home/jmendietaes/data/2021/chip/allProcessed/bamfiles/valid \
#/home/jmendietaes/data/2021/chip/allProcessed/furtherAnalysis/chipCorrelation


# script to compute correlation of BAMs for the same ChIP
bamsPath=$1
#bamsPath="/home/jmendietaes/data/2021/chip/allProcessed/bamfiles/valid"
outMatrixP=$2
#outMatrixP="/home/jmendietaes/data/2021/chip/allProcessed/furtherAnalysis/chipCorrelation"

# generate specific temporal folder for this job
tempFolder="${outMatrixP}/tempSubs_$((1 + $RANDOM % 1000000))"

# wether we want to compute the metrics adding all the control files in the comparison
CompareWithControl="NO"  # "YES" or "NO"

# load modules
module load SAMtools/1.12-GCC-10.2.0
module load Sambamba/0.7.0
export PATH="/home/jmendietaes/programas/miniconda3/bin:$PATH"
export PYTHONPATH=/home/jmendietaes/programas/miniconda3/bin/python3.8


##################### FUNCTIONS ####################
getMinReadNumber () {
    bamList=$1
    # lengths has to be local just in case
    local lengths=()
    # get reads in each of the samples
    for ba in ${bamList}; do
        count=$(samtools idxstats ${ba} | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {print total}');
        lengths[${#lengths[@]}]=${count}
    done
    # Get the minimum number of reads
    oldIFS=$IFS
    IFS=$'\n'
    minVal=$(echo "${lengths[*]}" | sort -n | head -1)
    IFS=$oldIFS
}

subSampleBams () {
    bamList=$1
    # store number of elements that were subsampled already
    prevLen="${#lengths[@]}"
    lengths=()
    # get reads in each of the samples
    for ba in ${bamList}; do
        count=$(samtools idxstats ${ba} | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {print total}');
        lengths[${#lengths[@]}]=${count}
    done
    # Get the minimum number of reads
    oldIFS=$IFS
    IFS=$'\n'
    minVal=$(echo "${lengths[*]}" | sort -n | head -1)
    IFS=$oldIFS
    echo -e "\n BAM file lengths:"
    echo ${bamList}
    echo "${lengths[*]}"

    # Get the number of reads to subsample
    factors=()
    for ba in ${bamList}; do
        FACTOR=$(samtools idxstats $ba | cut -f3 | awk -v COUNT=$minVal 'BEGIN {total=0} {total += $1} END {print COUNT/total}');
        factors[${#factors[@]}]=${FACTOR}
    done

    # check if the new minimum is the same as the one before
    if [[ ${minValPre} != ${minVal} ]]; then
        echo "Subsampling all files files"
        # clear previous files
        if [[ -e ${tempFolder}"/temp2_*"${Ig_prot}"*bam*" ]] ; then
            rm ${tempFolder}"/temp2_*"${Ig_prot}"*bam*"
        fi
        # They are different, so we need to subsample all again
        # subsample
        i=0
        subProtBams=()
        for ba in ${bamList}; do
            outba=`basename $ba`
            outba=${tempFolder}"/temp2_"${outba}
            sambamba view -h -t $SLURM_CPUS_PER_TASK -s ${factors[$i]} -f bam --subsampling-seed=12345 $ba -o $outba
            i=$(($i+1))
            subProtBams[${#subProtBams[@]}]=${outba}
        done
        subProtBams=`echo ${subProtBams[*]}`
    else
        # They are equal, so we only need to subsample the new files
        echo "Recycling previous files"
        # subsample
        i=0
        subProtBams=()
        for ba in ${bamList}; do
            # we subsample new files
            if [[ $i -ge ${prevLen} ]]; then
                outba=`basename $ba`
                outba=${tempFolder}"/temp2_"${outba}
                sambamba view -h -t $SLURM_CPUS_PER_TASK -s ${factors[$i]} -f bam --subsampling-seed=12345 $ba -o $outba
                i=$(($i+1))
                subProtBams[${#subProtBams[@]}]=${outba}
            # we recycle the previous subsampling
            else
                outba=`basename $ba`
                outba=${tempFolder}"/temp2_"${outba}
                i=$(($i+1))
                subProtBams[${#subProtBams[@]}]=${outba}
            fi
        done
        subProtBams=`echo ${subProtBams[*]}`
    fi
}

######################## MAIN CODE ##########################

if [ ! -e ${tempFolder} ]; then
    mkdir -p ${tempFolder}
fi


PROTS=($(cat $bamsPath/chipMarcs.txt))
Ig_prot=${PROTS[$SLURM_ARRAY_TASK_ID - 1]}

# get space separated paths of selected files
protBams=$(ls ${bamsPath}/*_${Ig_prot}*bam | tr '\n' ' ')

# set PCA symbols ‘<’,’>’,’o’)
symbols=''
for p in ${protBams}; do  
    pp=`basename ${p}`; 
    if [[ ${pp} == Mye* || ${pp} == ChIP* ]]; then 
        symbols="${symbols} <" 
    elif [[ ${pp} == DM*  ]]; then
        symbols="${symbols} >" 
    elif [[ ${pp} == GMPvitro*  ]]; then
        symbols="${symbols} s" 
    elif [[ ${pp} == MEPvitro*  ]]; then
        symbols="${symbols} p" 
    elif [[ ${pp} == LSKvitro*  ]]; then
        symbols="${symbols} D" 
    else
        symbols="${symbols} o" 
    fi 
done


## compare bams without controls
labels=`for i in $protBams; do basename ${i} | cut -d '.' -f 1 | cut -d '_' -f 1,2,3; done | tr '\n' ' '`

# check if this comparison was already done (leave plots in case there were code changes)
doneComp=`head -2 ${outMatrixP}/${Ig_prot}_PCAinfo.txt | tail -1 | awk '{$1=""; $(NF)=""; print $0}'`
# Define needed variables
lengths=()
minValPre='No'
# get minimum read number from our sample files anc chek if files already exist
getMinReadNumber "${protBams}"
matrixOut="${outMatrixP}/${Ig_prot}_${minVal}Read_coverageComp.npz"

if [[ ! (-e ${matrixOut}) ||  `echo $labels` != `echo $doneComp` ]]; then 
    ## To do a fair comparison we need to subsample bams and leave them with similar number of reads
    # call subsampling function
    subSampleBams "${protBams}"

    ## Then we compute the matrix
    multiBamSummary bins --bamfiles ${subProtBams} -o ${matrixOut} \
                    --labels ${labels} --numberOfProcessors ${SLURM_CPUS_PER_TASK}
fi

plotCorrelation --corData ${matrixOut} --corMethod spearman --whatToPlot heatmap \
            -o ${outMatrixP}/${Ig_prot}_${minVal}Read_coverageComp_heatmap.pdf --skipZeros --plotTitle ${Ig_prot} \
            --plotNumbers --colorMap viridis

plotPCA --corData ${matrixOut} -o ${outMatrixP}/${Ig_prot}_${minVal}Read_coverageComp_PCA.pdf \
    --plotTitle ${Ig_prot} --outFileNameData ${outMatrixP}/${Ig_prot}_${minVal}Read_PCAinfo.txt \
    --markers ${symbols}



## compare bams with controls
if [[ $CompareWithControl == 'YES' ]]; then

    # only run in non-control samples
    if [[ $Ig_prot != 'IgG' && $Ig_prot != 'input' ]]; then 
        # get IgG controls
        iggBams=$(ls ${bamsPath}/*_IgG*bam | tr '\n' ' ')
        # get inputs
        inputBams=$(ls ${bamsPath}/*_input*bam | tr '\n' ' ')
        # concatenate
        compareBams="${protBams} ${iggBams} ${inputBams}"
    elif [[ $Ig_prot != 'IgG' ]]; then
        # get inputs
        inputBams=$(ls ${bamsPath}/*_input*bam | tr '\n' ' ')
        # concatenate
        compareBams="${protBams} ${inputBams}"
    elif [[ $Ig_prot != 'input' ]]; then
        # get IgG controls
        iggBams=$(ls ${bamsPath}/*_IgG*bam | tr '\n' ' ')
        # concatenate
        compareBams="${protBams} ${iggBams}"
    fi

    # prepare additional parameters
    labels=`for i in $compareBams; do basename ${i} | cut -d '.' -f 1 | cut -d '_' -f 1,2,3; done | tr '\n' ' '`

    # check if this comparison was already done (ley plots in case there were code changes)
    doneComp=`head -2 ${outMatrixP}/${Ig_prot}_PCAinfo_wcontrol.txt | tail -1 | awk '{$1=""; $(NF)=""; print $0}' | tr '\n' ' '`

    # get minimum read number from our sample files anc chek if files already exist
    minValPre=$minVal
    getMinReadNumber "${compareBams}"
    matrixOut="${outMatrixP}/${Ig_prot}_${minVal}Read_coverageComp_wcontrol.npz"
    if [[ ! (-e ${matrixOut}) ||  `echo $labels` != `echo $doneComp` ]]; then 
        ## subsample bams
        subSampleBams "${compareBams}"

        ## Then we compute the matrix
        multiBamSummary bins --bamfiles ${subProtBams} -o ${matrixOut} \
                        --labels ${labels} --numberOfProcessors ${SLURM_CPUS_PER_TASK}
    fi

    

    plotCorrelation --corData ${matrixOut} --corMethod spearman --whatToPlot heatmap \
                -o ${outMatrixP}/${Ig_prot}_${minVal}Read_coverageComp_wcontrol_heatmap.pdf --skipZeros --plotTitle ${Ig_prot} \
                --plotNumbers --colorMap viridis

    # add PCA symbol for controls
    symbolsC=$symbols
    ninput=$((`echo $iggBams | wc -w`+`echo $inputBams | wc -w`))
    for s in `seq 1 ${ninput}`; do
        symbolsC="${symbolsC} o"
    done

    plotPCA --corData ${matrixOut} -o ${outMatrixP}/${Ig_prot}_${minVal}Read_coverageComp_wcontrol_PCA.pdf \
        --plotTitle ${Ig_prot} --outFileNameData ${outMatrixP}/${Ig_prot}_${minVal}Read_PCAinfo_wcontrol.txt \
        --markers ${symbolsC}
fi

# clear temporal files
rm -r ${tempFolder}