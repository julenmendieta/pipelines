#!/bin/bash
# -*- ENCODING: UTF-8 -*-

##===============================================================================
## SLURM VARIABLES
#SBATCH --job-name=peakAnalysis
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=1-01:00:00
#SBATCH -p medium
#SBATCH -o /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.out  
#SBATCH -e /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.err 


# HOW TO RUN ME
#sbatch /home/jmendietaes/programas/PhD/ChIP/cluster/02b_peakAnalysis_I.sh \
#/home/jmendietaes/data/2021/chip/allProcessed \


# NOTES
# Peak caling of merges
# CREO QUE A PARTIR DE AQUI VOY A HACER OTRO SCRIPT, PARA ELLO TENDRE QUE GUARDAR EL BAMFILE 
# EN ALGUN SITIO EN CONCRETO DONDE COINCIDAN LAS REPLICAS DE LA MISMA MUESTRA
# pero tambien voy a sacar los picos de las muestras por separado
# a MACS2 se le puede pasar los BAM en separado y el los junta


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
# species shortcut for MACS
species="mm"
speciesGenome="mm10"

# extend variables
bamsPath="${basePath}/bamfiles/valid"
outpath=${basePath}"/furtherAnalysis"


allbams=$(find ${bamsPath}/*bam -printf "${bamsPath}/%f ")
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

#############################
#PEAK CALLING: MACS2
#############################


echo -e "Starting Peak calling -----------------------------------------------\n"


# Create output dir
if [ ! -e ${outpath}/peakCalling/MACS2/logs ]; then
    mkdir -p ${outpath}/peakCalling/MACS2/logs
    mkdir -p ${outpath}/peakCalling/MACS2/peaks/summary
fi

cd ${outpath}/peakCalling/MACS2/peaks

# peaktype="narrowPeak"
#peaktype="broadPeak"
# if [ "$peaktype" == "broad" ]; then
#     peaktype2="--broad"
# else
#     peaktype2=''
# fi

## we go for each cell
for bam in ${allbams}; do
    cell=$(basename $bam | cut -d '.' -f 1 | cut -d '_' -f 1)
    cellControls=`echo ${allbams} | tr ' ' '\n' | grep "${cell}_" | \
                    { grep -e "input" -e "IgG" || :; }`
    label=$(basename $bam | cut -d '.' -f 1 | cut -d '_' -f 1,2,3)
    # for now we use the first control to normalise signal (by sort should be IgG)
    controlbam=`echo ${cellControls} | cut -d ' ' -f 1`
    summaryFile="${outpath}/peakCalling/MACS2/peaks/summary/summary_${label}.txt"

    total_reads="empty"

    # narrow peaks
    peaktype='narrowPeak'

    # check if the file exists of it was created with a previous bam version 
    fileNotExistOrOlder "${outpath}/peakCalling/MACS2/peaks/${label}_peaks_${peaktype}.xls" "${bam} ${controlbam}"
    # this outputs analyse as yes or no in lowercase

    if [[ ${analyse} == "yes" ]]; then
        total_reads=$(samtools view -c ${bam})

        macs2 callpeak \
                -t ${bam} \
                -c ${controlbam} \
                -f BAMPE \
                -g $species \
                -n $label \
                --keep-dup all \
                --outdir ${outpath}/peakCalling/MACS2/peaks 2> ${outpath}/peakCalling/MACS2/logs/${label}_macs2_${peaktype}.log

        mv ${outpath}/peakCalling/MACS2/peaks/${label}_peaks.xls ${outpath}/peakCalling/MACS2/peaks/${label}_peaks_${peaktype}.xls

        npeaks=$(cat ${outpath}/peakCalling/MACS2/peaks/${label}_peaks.${peaktype} | wc -l)
        reads_in_peaks=$(bedtools sort -i ${outpath}/peakCalling/MACS2/peaks/${label}_peaks.${peaktype} \
            | bedtools merge -i stdin | bedtools intersect -u -nonamecheck \
            -a ${bam} -b stdin -ubam | samtools view -c)
        FRiP=$(awk "BEGIN {print "${reads_in_peaks}"/"${total_reads}"}")
        # report
        echo -e "NUMBER OF NARROW PEAKS \t ${npeaks} \n" >> ${summaryFile}
        echo -e "total_reads \t reads_in_peaks \t FRIP \n" >> ${summaryFile}
        echo -e "${total_reads} \t ${reads_in_peaks} \t ${FRiP}" >> ${summaryFile}
    fi

    # broad peaks
    peaktype="broadPeak"

    # check if the file exists of it was created with a previous bam version 
    fileNotExistOrOlder "${outpath}/peakCalling/MACS2/peaks/${label}_peaks_${peaktype}.xls" "${bam} ${controlbam}"
    # this outputs analyse as yes or no in lowercase

    if [[ ${analyse} == "yes" ]]; then
        if [[ $total_reads == "empty" ]]; then
            total_reads=$(samtools view -c ${bam})
        fi

        macs2 callpeak \
                -t ${bam} \
                -c ${controlbam} \
                --broad \
                -f BAMPE \
                -g $species \
                -n $label \
                --keep-dup all \
                --outdir ${outpath}/peakCalling/MACS2/peaks 2> ${outpath}/peakCalling/MACS2/logs/${label}_macs2_${peaktype}.log
        
        # -q 0.05 as default (at least for broad)
        mv ${outpath}/peakCalling/MACS2/peaks/${label}_peaks.xls ${outpath}/peakCalling/MACS2/peaks/${label}_peaks_${peaktype}.xls

        npeaks=$(cat ${outpath}/peakCalling/MACS2/peaks/${label}_peaks.${peaktype} | wc -l)
        reads_in_peaks=$(bedtools sort -i ${outpath}/peakCalling/MACS2/peaks/${label}_peaks.${peaktype} \
            | bedtools merge -i stdin | bedtools intersect -u -nonamecheck \
            -a ${bam} -b stdin -ubam | samtools view -c)
        FRiP=$(awk "BEGIN {print "${reads_in_peaks}"/"${total_reads}"}")
        # report
        echo -e "NUMBER OF BROAD PEAKS \t ${npeaks} \n" >> ${summaryFile}
        echo -e "total_reads \t reads_in_peaks \t FRIP \n" >> ${summaryFile}
        echo -e "${total_reads} \t ${reads_in_peaks} \t ${FRiP}" >> ${summaryFile}
    fi
    # wget https://raw.githubusercontent.com/nf-core/chipseq/master/assets/multiqc/peak_count_header.txt
    # wget https://raw.githubusercontent.com/nf-core/chipseq/master/assets/multiqc/frip_score_header.txt
    # peak_count_header="${extraFilePath}/peak_count_header.txt"
    # frip_score_header="${extraFilePath}/frip_score_header.txt"

    # ipflagstat=/home/jmendietaes/data/2021/chip/allProcessed/furtherAnalysis/MACS2/flagstat.txt

    # echo $npeaks | awk -v ip=$ip -v OFS='\t' '{ print ip, $1 }' | cat $peak_count_header - > ${ip}_${peaktype}s.count_mqc.tsv

    # Maybe im too flexible, i could add -f 0.10 or -F 0.10 to ask at least for 10% overlap
    # awk -v FRiP=$FRiP -v OFS='\t' '{print ip, FRiP}' | cat $frip_score_header - > ${ip}_${peaktype}s.FRiP_mqc.tsv
done

echo -e "Peak calling - finished -----------------------------------------------\n"


###############################
# HOMER: annotate peaks
##############################

# anotate peaks HOMER
# http://homer.ucsd.edu/homer/ngs/annotation.html

# Create output dir
if [ ! -e ${outpath}/HOMER/peakAnnotation ]; then
    mkdir -p ${outpath}/HOMER/peakAnnotation
fi

echo -e "Starting HOMER annotation -----------------------------------------------\n"


for bam in ${allbams}; do
    for peaktype in 'narrowPeak' 'broadPeak'; do

        label=$(basename $bam | cut -d '.' -f 1 | cut -d '_' -f 1,2,3)
        annotationOut=${outpath}/HOMER/peakAnnotation/${label}_${peaktype}.annotatePeaks.txt
        peakIn=${outpath}/peakCalling/MACS2/peaks/${label}_peaks.${peaktype}
        
        # check if the file exists or it was created with a previous peaks version 
        fileNotExistOrOlder "${annotationOut}" "${peakIn}"
        # this outputs analyse as yes or no in lowercase
        if [[ ${analyse} == "yes" ]]; then
            annotatePeaks.pl \
                    ${peakIn} \
                    ${speciesGenome} \
                    -gid \
                    -cpu ${SLURM_CPUS_PER_TASK} \
                    -annStats ${outpath}/HOMER/peakAnnotation/${label}_${peaktype}.annotateStats.txt
                    > ${annotationOut}
        fi
    done
done

echo -e "HOMER annotation - Finished -----------------------------------------------\n"

###############################
# QC: MACS2 quality check plots
##############################
# plot metrics
if [ ! -e ${outpath}/peakCalling/MACS2/QC ]; then
    mkdir -p ${outpath}/peakCalling/MACS2/QC
fi

echo -e "Starting peak QC -----------------------------------------------\n"


allCell=`for filename in ${allLabels}; do 
            mapLib=(${filename//_/ }); 
            echo ${mapLib[0]}; done | sort | uniq`

for cell in ${allCell}; do
    for peaktype in 'narrowPeak' 'broadPeak'; do
        echo "${cell}_${peaktype}_macs2"
        peaksPath=${outpath}/peakCalling/MACS2/peaks
        # get comma separated string with all the peak files to compare
        peakfiles=$(find ${peaksPath}/${cell}_*peaks.${peaktype} -printf "${peaksPath}/%f,")
        labels=$(find ${peaksPath}/${cell}_*peaks.${peaktype} -printf "%f " -exec basename {} \; )
        labels=$(for la in $labels; do echo $la | sed "s/_peaks.${peaktype}//g"; done | tr '\n' ',' )

        # check if the file exists or it was created with a previous peaks version 
        peakfilesSpace=$(find ${peaksPath}/${cell}_*peaks.${peaktype} -printf "${peaksPath}/%f ")
        fileNotExistOrOlder "${outpath}/peakCalling/MACS2/QC/${cell}_${peaktype}_macs2.plots.pdf" "${peakfilesSpace}"
        # this outputs analyse as yes or no in lowercase
        if [[ ${analyse} == "yes" ]]; then
            # https://github.com/nf-core/chipseq
            Rscript ${scriptsPath}/ChIP/cluster/02_NR_plot_macs_qc.R \
                    -i ${peakfiles} \
                    -s ${labels} \
                    -o ${outpath}/peakCalling/MACS2/QC \
                    -p "${cell}_${peaktype}_macs2" # out prefix
        fi

        # get comma separated string with all the peak annotation files to compare (avoid controls)
        annotPath=${outpath}/HOMER/peakAnnotation
        annotfiles=$(find ${annotPath}/${cell}_*${peaktype}.annotatePeaks.txt | \
                    { grep -v -e "input" -v -e "IgG" || :; } | tr '\n' ',' )
        labels=$(find ${annotPath}/${cell}_*${peaktype}.annotatePeaks.txt -exec basename {} \; | \
                    { grep -v -e "input" -v -e "IgG" || :; } )
        labels=$(for la in $labels; do echo $la | sed "s/_${peaktype}.annotatePeaks.txt//g"; done | tr '\n' ',' )

        # check if the file exists or it was created with a previous peaks version 
        annotfilesSpace=$(find ${annotPath}/${cell}_*${peaktype}.annotatePeaks.txt -printf "${annotPath}/%f ")
        fileNotExistOrOlder "${outpath}/peakCalling/MACS2/QC/${cell}_${peaktype}_macs2_annotatePeaks.plots.pdf" "${annotfilesSpace}"
        # this outputs analyse as yes or no in lowercase
        if [[ ${analyse} == "yes" ]]; then
            Rscript ${scriptsPath}/ChIP/cluster/02_NR_plot_homer_annotatepeaks.R \
                    -i ${annotfiles} \
                    -s ${labels} \
                    -o ${outpath}/peakCalling/MACS2/QC \
                    -p "${cell}_${peaktype}_macs2_annotatePeaks"
        fi
    done
done

echo -e "peak QC - Finished -----------------------------------------------\n"


#########################
# CONSENSUS PEAKS ANALYSIS: Consensus peaks across samples, create boolean 
#     filtering file, SAF file for featureCounts and UpSetR plot for intersection
#########################

## We first do it in the whole dataset as a trial, but in the future ill do it only between the biological replicates
if [ ! -e ${outpath}/peakCalling/MACS2/consensusPeaks ]; then
	mkdir -p ${outpath}/peakCalling/MACS2/consensusPeaks
fi

echo -e "Starting consensus peak analysis -------------------------------------\n"


chip="allmerged"
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
                { grep -v -e "input" -v -e "IgG" || :; } )
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



## Same by comparing samples of same chip
echo -e "Starting same-chip peak comparison -------------------------------------\n"

sameChipCons=${outpath}/peakCalling/MACS2/consensusPeaks/bySameChip
if [ ! -e ${sameChipCons} ]; then
	mkdir -p ${sameChipCons}
fi

for peaktype in narrowPeak broadPeak; do
    # select the columns to peak in each peak calling case
    if [ ${peaktype} == "narrowPeak" ]; then
        mergecols=`seq 2 10 | tr '\n' ','`
        expandparam='--is_narrow_peak'
    elif [ ${peaktype} == "broadPeak" ]; then
        mergecols=`seq 2 9 | tr '\n' ','`
        expandparam=''
    fi

    # get list with all checked chip
    allChip=$(\
    for filename in ${bamsPath}/*bam; do 
        mapLib=(${filename//_/ }); 
        mapLib=${mapLib[1]}; 
        mapLib=(${mapLib//-/ }); 
        echo ${mapLib[0]}; done | \
        grep -v input | grep -v IgG | sort | uniq)

    
    # compare the ones for which we either have replicates or other conditions
    for chip in ${allChip}; do
        # here we check if we have more than one file for that chip
        peakFiles=$(find -L ${outpath}/peakCalling/MACS2/peaks/*.${peaktype} -maxdepth 1  -type f ! -size 0 | \
                { grep -e ${chip} || :; } )
        fileLabels=$(for f in $peakFiles; do echo ${f##*/} |sed "s/_peaks.${peaktype}//g"; done | tr '\n' ',')
        prefix="${chip}_${peaktype}_consensusPeaks"


        nfiles=$(echo $peakFiles | wc -w)
        if [[ ${nfiles} -gt 1 ]]; then    

            RoutPlot=${sameChipCons}/${prefix}.boolean.intersect.pdf
            fileNotExistOrOlder "${RoutPlot}" "${peakFiles}"
            # this outputs analyse as yes or no in lowercase
            if [[ ${analyse} == "yes" ]]; then

                sort -T '.' -k1,1 -k2,2n ${peakFiles} \
                | mergeBed -c $mergecols -o collapse > ${sameChipCons}/${prefix}.txt

                python ${scriptsPath}/ChIP/cluster/02_NR_macs2_merged_expand.py ${sameChipCons}/${prefix}.txt \
                    ${fileLabels} \
                    ${sameChipCons}/${prefix}.boolean.txt \
                    $expandparam

                consensusPeakBed=${sameChipCons}/${prefix}.bed
                awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print $1, $2, $3, $4, "0", "+" }' \
                    ${sameChipCons}/${prefix}.boolean.txt > ${consensusPeakBed}
                echo -e "GeneID\tChr\tStart\tEnd\tStrand" > ${sameChipCons}/${prefix}.saf
                awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print $4, $1, $2, $3,  "+" }' \
                    ${sameChipCons}/${prefix}.boolean.txt >> ${sameChipCons}/${prefix}.saf
                Rscript ${scriptsPath}/ChIP/cluster/02_NR_plot_peak_intersect.r \
                    -i ${sameChipCons}/${prefix}.boolean.intersect.txt \
                    -o ${RoutPlot}
            fi
        
        fi

    done

done


echo -e "same-chip peak comparison - Finished ---------------------------\n"



############################################################
# Annotate consensus peaks with HOMER, and add annotation to boolean output file
############################################################

if [ ! -e ${outpath}/HOMER/consensusPeaks ]; then
	mkdir -p ${outpath}/HOMER/consensusPeaks
fi

echo -e "Starting consensus peak annotations ------------------------------\n"


chip="allmerged"
for peaktype in narrowPeak broadPeak; do

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


## Same by comparing samples of same chip

echo -e "Starting consensus same-chip peak annotations -----------------------\n"

sameChipHomer=${outpath}/HOMER/consensusPeaks/bySameChip
if [ ! -e ${sameChipHomer} ]; then
	mkdir -p ${sameChipHomer}
fi


# get all analysed ChIP
more1Chip=$(\
    for filename in ${sameChipCons}/*bed; do 
        filename=$(basename ${filename})
        mapLib=(${filename//_/ }); 
        echo ${mapLib[0]}; done | \
        sort | uniq)

for chip in ${more1Chip}; do
    for peaktype in narrowPeak broadPeak; do

        prefix="${chip}_${peaktype}_consensusPeaks"
        consensusPeakBed=${sameChipCons}/${prefix}.bed

        # check if the file exists or it was created with a previous peaks version 
        boolAnotMatr=${sameChipHomer}/${prefix}.boolean.annotatePeaks.txt
        fileNotExistOrOlder "${boolAnotMatr}" "${consensusPeakBed}"
        # this outputs analyse as yes or no in lowercase
        if [[ ${analyse} == "yes" ]]; then

            annotatePeaks.pl \
                    ${consensusPeakBed} \
                    ${speciesGenome} \
                    -gid \
                    -cpu ${SLURM_CPUS_PER_TASK} \
                    -annStats ${sameChipHomer}/${prefix}.annotateStats.txt \
                    > ${sameChipHomer}/${prefix}.annotatePeaks.txt

            cut -f2- ${sameChipHomer}/${prefix}.annotatePeaks.txt | \
                awk 'NR==1; NR > 1 {print $0 | "sort -T '.' -k1,1 -k2,2n"}' | \
                cut -f6- > ${sameChipHomer}/tmp.txt
            paste ${sameChipCons}/${prefix}.boolean.txt \
                ${sameChipHomer}/tmp.txt > ${boolAnotMatr}
        fi
    done
done

echo -e "consensus same-chip peak annotations - Finished ---------------\n"


###########################################################
# Count reads in consensus peaks with featureCounts
##########################################################

featureCpath=${outpath}/peakCalling/MACS2/consensusPeaks/featureCounts
if [ ! -e ${featureCpath} ]; then
	mkdir -p ${featureCpath}
fi
cd {featureCpath}


chip="allmerged"
for peaktype in narrowPeak broadPeak; do
    prefix="${chip}_${peaktype}_consensusPeaks"
    peakFiles=$(find -L ${outpath}/peakCalling/MACS2/peaks/*.${peaktype} -maxdepth 1  -type f ! -size 0 | \
                { grep -v -e "input" -v -e "IgG" || :; } )
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


## Same by comparing samples of same chip

featureCPairpath=${outpath}/peakCalling/MACS2/consensusPeaks/bySameChip/featureCounts
if [ ! -e ${featureCPairpath} ]; then
	mkdir -p ${featureCPairpath}
fi
cd ${featureCPairpath}


for chip in ${more1Chip}; do
    for peaktype in narrowPeak broadPeak; do
        prefix="${chip}_${peaktype}_consensusPeaks"
        peakFiles=$(find -L ${outpath}/peakCalling/MACS2/peaks/*.${peaktype} -maxdepth 1  -type f ! -size 0 | \
                    { grep "${chip}" || :; } )
        fileLabels=$(for f in $peakFiles; do echo ${f##*/} |sed "s/_peaks.${peaktype}//g"; done)
        # get a list with the bam files in the same order as peak files
        bamfiles=$(for f in ${fileLabels}; do find -L ${bamsPath}/${f}*bam -printf "${bamsPath}/%f "; 
                        done | tr '\n' ' ')

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

###########################################################
# Read count to CPM
###########################################################

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
                    -S "allProcessed.bamfiles.valid." \
                    -E ".sort.rmdup.rmblackls.rmchr.bam" 

    fi
done

## Same by comparing samples of same chip

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
            sname=`basename $bam  | sed 's/.sort.rmdup.rmblackls.rmchr.bam//g'`
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
                        -d ${prevFields} \
                        -S "allProcessed.bamfiles.valid." \
                        -E ".sort.rmdup.rmblackls.rmchr.bam" 

        fi
    done
done



###########################################################
# Differential analysis with DESeq2
###########################################################

# if [ ! -e ${outpath}/DESeq2/MACS2 ]; then
# 	mkdir -p ${outpath}/DESeq2/MACS2
# fi

# chip="allmerged"
# for peaktype in narrowPeak broadPeak; do
#     prefix="${chip}_${peaktype}_consensusPeaks" 

#     featureOut=${featureCpath}/${prefix}.featureCounts.txt


#     Rscript ${scriptsPath}/ChIP/cluster/02_NR_featurecounts_deseq2.r \
#             --featurecount_file ${featureOut} \
#             --bam_suffix '.sort.rmdup.rmblackls.rmchr.bam' \
#             --outdir ${outpath}/DESeq2/ \
#             --outprefix $prefix \
#             --outsuffix '' \
#             --cores ${SLURM_CPUS_PER_TASK}

#     sed 's/deseq2_pca/deseq2_pca_${task.index}/g' <$deseq2_pca_header >tmp.txt
#     sed -i -e 's/DESeq2 /${antibody} DESeq2 /g' tmp.txt
#     cat tmp.txt ${prefix}.pca.vals.txt > ${prefix}.pca.vals_mqc.tsv
#     sed 's/deseq2_clustering/deseq2_clustering_${task.index}/g' <$deseq2_clustering_header >tmp.txt
#     sed -i -e 's/DESeq2 /${antibody} DESeq2 /g' tmp.txt
#     cat tmp.txt ${prefix}.sample.dists.txt > ${prefix}.sample.dists_mqc.tsv
# done


# ## Same by comparing samples of same chip
# if [ ! -e ${outpath}/DESeq2/MACS2/bySameChip ]; then
# 	mkdir -p ${outpath}/DESeq2/MACS2/bySameChip
# fi

# for chip in ${more1Chip}; do
#     for peaktype in narrowPeak broadPeak; do
#         prefix="${chip}_${peaktype}_consensusPeaks" 

#         featureOut=${featureCpath}/${prefix}.featureCounts.txt


#         Rscript ${scriptsPath}/ChIP/cluster/02_NR_featurecounts_deseq2.r \
#                 --featurecount_file ${featureOut} \
#                 --bam_suffix '.sort.rmdup.rmblackls.rmchr.bam' \
#                 --outdir ${outpath}/DESeq2/ \
#                 --outprefix $prefix \
#                 --outsuffix '' \
#                 --cores ${SLURM_CPUS_PER_TASK}

#         sed 's/deseq2_pca/deseq2_pca_${task.index}/g' <$deseq2_pca_header >tmp.txt
#         sed -i -e 's/DESeq2 /${antibody} DESeq2 /g' tmp.txt
#         cat tmp.txt ${prefix}.pca.vals.txt > ${prefix}.pca.vals_mqc.tsv
#         sed 's/deseq2_clustering/deseq2_clustering_${task.index}/g' <$deseq2_clustering_header >tmp.txt
#         sed -i -e 's/DESeq2 /${antibody} DESeq2 /g' tmp.txt
#         cat tmp.txt ${prefix}.sample.dists.txt > ${prefix}.sample.dists_mqc.tsv
#     done
# done

###########################################################
# IGV
###########################################################


exit 0