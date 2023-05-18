#!/bin/bash
# -*- ENCODING: UTF-8 -*-

##===============================================================================
## SLURM VARIABLES
#SBATCH --job-name=ATACanalysis
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=05:00:00
#SBATCH -p short
#SBATCH -o /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.out  
#SBATCH -e /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.err 


# HOW TO RUN ME
#sbatch /home/jmendietaes/programas/pipelines/ATAC/cluster/02b_peakAnalysis_I.sh \
#/home/jmendietaes/data/2021/ATAC/allProcessed

# OBJECTIVE
# call peaks, annotate, make consensus peaks, get reads in peaks, CPM and merge
# all in same table
# This script has been recycled from ChIP-seq peak calling pipeline and retains
# some unused checks

# path where we have the folder structure for chip analysis 
# (inside we have ${basePath}/bamfiles/valid/)
basePath=$1
#basePath="/home/jmendietaes/data/2021/ATAC/allProcessed"
# Location of the analysis header files downloaded from 
# https://github.com/nf-core/chipseq/tree/master/assets/multiqc
#extraFilePath=$2
#extraFilePath="/home/jmendietaes/data/2021/chip/analysisFiles"

# Where to look for bam files and where to store output tree
bamsPath="${basePath}/bamfiles/valid/08_paperChip"
outpath=${basePath}"/furtherAnalysis/08_paperChip"


# FILE NAMING FORMAT
# [cellType]_[chip]_[date]_[extra?].[bamfilteringKeys].bam
# for final table naming only first 3 sections between '_'
#   will be used, but only 2 are also accepted
# chip can contain sections separated by '-', and only the
#   first section will be used to count ChIPs as same
# extra sections wont be used in final naming
# bamfilteringKeys are supposed to be separated by dots

# GTF file for annotation (top be consistent with scRNA data)
# Set to FALSE if you wnat HOMER's default UCSC refGene annotation
#gtfFile=/home/jmendietaes/data/2021/singleCell/additionalFiles/refdata-gex-mm10-2020-A/genes/genes.gtf
gtfFile=FALSE

# If we want a column focussed on repeated elements only
# Path to Homer file with repeat element locations
repeatsPath=/beegfs/easybuild/CentOS/7.5.1804/Skylake/software/Homer/4.10-foss-2018b/data/genomes/mm10/mm10.repeats
#repeatsPath=FALSE

# path for the location of the pipeline scripts
scriptsPath="/home/jmendietaes/programas/pipelines"
# species shortcut for MACS
species="mm"
speciesGenome="mm10"
# Path to used genomic reference
REFERENCE_DIR="/home/jmendietaes/referenceGenomes/mm10_reordered/mm10.reordered"
chr_genome_size=$REFERENCE_DIR".sizes"

# Path to bedgraphToBigwig script
bedGraphToBigWig="/home/jmendietaes/programas/pipelines/general/bedGraphToBigWig"

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

# get extra parameters for annotation
if [[ $gtfFile == "FALSE" ]]; then
    extraAnnot=""
else
    extraAnnot="-gtf ${gtfFile}"
fi

#############################
#PEAK CALLING: MACS2
#############################


echo -e "Starting Peak calling -----------------------------------------------\n"


# Create output dir
if [ ! -e ${outpath}/peakCalling/MACS2/logs ]; then
    mkdir -p ${outpath}/peakCalling/MACS2/logs
    mkdir -p ${outpath}/peakCalling/MACS2/bigwig
    mkdir -p ${outpath}/peakCalling/MACS2/peaks/summary
    mkdir -p ${outpath}/peakCalling/MACS2/beds
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
    # proceed to call peaks for no controls
    if ! { [[ $bam == *"_IgG"* ]] || [[ $bam == *"_input"*  ]]; } ; then
        
        cell=$(basename $bam | cut -d '.' -f 1 | cut -d '_' -f 1)

        label=$(basename $bam | cut -d '.' -f 1 | cut -d '_' -f 1,2,3)

        summaryFile="${outpath}/peakCalling/MACS2/peaks/summary/summary_${label}.txt"

        total_reads="empty"

        # echo file and control
        echo "Bams used for peak calling:"
        echo $bam
        echo

        # narrow peaks
        peaktype='narrowPeak'

        # check if the file exists of it was created with a previous bam version 
        fileNotExistOrOlder "${outpath}/peakCalling/MACS2/peaks/${label}_peaks_${peaktype}.xls" "${bam}"
        # this outputs analyse as yes or no in lowercase

        outbed=$(basename $bam)
        outbed="${outpath}/peakCalling/MACS2/beds/${outbed::-4}.bed"

        if [[ ${analyse} == "yes" ]]; then
            total_reads=$(samtools view -c ${bam})

            # We convert the BAM file to BED format because when we set the 
            #extension size in MACS2, it will only consider one read of the 
            #pair while here we would like to use the information from both
            # # Convert bam to bed
            bamToBed -i ${bam} > ${outbed}

            # We want to focuss in the cut sites of Tn5, butIf we only assess 
            #the coverage of the 5’ extremity of the reads, 
            #the data would be too sparse and it would be impossible to call 
            #peaks. Thus, we will extend the start sites of the reads by 150bp 
            #(75bp in each direction) to assess coverage.
            macs2 callpeak \
                    -t ${outbed} \
                    -f BED \
                    -g $species \
                    -n $label \
                    --keep-dup all \
                    --nomodel --shift -75 --extsize 150 \
                    --outdir ${outpath}/peakCalling/MACS2/peaks 2> \
                        ${outpath}/peakCalling/MACS2/logs/${label}_macs2_${peaktype}.log

            mv ${outpath}/peakCalling/MACS2/peaks/${label}_peaks.xls ${outpath}/peakCalling/MACS2/peaks/${label}_peaks_${peaktype}.xls
            #mv ${outpath}/peakCalling/MACS2/peaks/${label}_treat_pileup.bdg \
            #   ${outpath}/peakCalling/MACS2/peaks/${label}_treat_pileup_${peaktype}.bdg

            npeaks=$(cat ${outpath}/peakCalling/MACS2/peaks/${label}_peaks.${peaktype} | wc -l)
            reads_in_peaks=$(bedtools sort -i ${outpath}/peakCalling/MACS2/peaks/${label}_peaks.${peaktype} \
                | bedtools merge -i stdin | bedtools intersect -u -nonamecheck \
                -a ${bam} -b stdin -ubam | samtools view -c)
            FRiP=$(awk "BEGIN {print "${reads_in_peaks}"/"${total_reads}"}")
            # report
            echo -e "NUMBER OF NARROW PEAKS \t ${npeaks} \n" >> ${summaryFile}
            echo -e "total_reads \t reads_in_peaks \t FRIP \n" >> ${summaryFile}
            echo -e "${total_reads} \t ${reads_in_peaks} \t ${FRiP}" >> ${summaryFile}

            # get bigwig
            #${bedGraphToBigWig} ${outpath}/peakCalling/MACS2/peaks/${label}_treat_pileup_${peaktype}.bdg $chr_genome_size \
            #                        ${outpath}/peakCalling/MACS2/bigwig/${label}_treat_pileup_${peaktype}.bw
            
        fi

        # broad peaks
        peaktype="broadPeak"

        # check if the file exists of it was created with a previous bam version 
        fileNotExistOrOlder "${outpath}/peakCalling/MACS2/peaks/${label}_peaks_${peaktype}.xls" "${bam}"
        # this outputs analyse as yes or no in lowercase

        if [[ ${analyse} == "yes" ]]; then
            if [[ $total_reads == "empty" ]]; then
                total_reads=$(samtools view -c ${bam})
            fi

            # If we only assess the coverage of the 5’ extremity of the reads, 
            #the data would be too sparse and it would be impossible to call 
            #peaks. Thus, we will extend the start sites of the reads by 150bp 
            #(75bp in each direction) to assess coverage.
            macs2 callpeak \
                    -t ${outbed} \
                    --broad \
                    -f BED \
                    -g $species \
                    -n $label \
                    --keep-dup all \
                    --nomodel --shift -75 --extsize 150 \
                    --outdir ${outpath}/peakCalling/MACS2/peaks \
                    2> ${outpath}/peakCalling/MACS2/logs/${label}_macs2_${peaktype}.log
            
            # -q 0.05 as default (at least for broad)
            mv ${outpath}/peakCalling/MACS2/peaks/${label}_peaks.xls \
                ${outpath}/peakCalling/MACS2/peaks/${label}_peaks_${peaktype}.xls
            #mv ${outpath}/peakCalling/MACS2/peaks/${label}_treat_pileup.bdg \
            #    ${outpath}/peakCalling/MACS2/peaks/${label}_treat_pileup_${peaktype}.bdg

            npeaks=$(cat ${outpath}/peakCalling/MACS2/peaks/${label}_peaks.${peaktype} | \
                    wc -l)
            reads_in_peaks=$(bedtools sort -i ${outpath}/peakCalling/MACS2/peaks/${label}_peaks.${peaktype} \
                | bedtools merge -i stdin | bedtools intersect -u -nonamecheck \
                -a ${bam} -b stdin -ubam | samtools view -c)
            FRiP=$(awk "BEGIN {print "${reads_in_peaks}"/"${total_reads}"}")
            # report
            echo -e "NUMBER OF BROAD PEAKS \t ${npeaks} \n" >> ${summaryFile}
            echo -e "total_reads \t reads_in_peaks \t FRIP \n" >> ${summaryFile}
            echo -e "${total_reads} \t ${reads_in_peaks} \t ${FRiP}" >> ${summaryFile}

            # get bigwig
            #${bedGraphToBigWig} ${outpath}/peakCalling/MACS2/peaks/${label}_treat_pileup_${peaktype}.bdg $chr_genome_size \
            #                        ${outpath}/peakCalling/MACS2/bigwig/${label}_treat_pileup_${peaktype}.bw

        fi
        # wget https://raw.githubusercontent.com/nf-core/chipseq/master/assets/multiqc/peak_count_header.txt
        # wget https://raw.githubusercontent.com/nf-core/chipseq/master/assets/multiqc/frip_score_header.txt
        # peak_count_header="${extraFilePath}/peak_count_header.txt"
        # frip_score_header="${extraFilePath}/frip_score_header.txt"

        # ipflagstat=/home/jmendietaes/data/2021/chip/allProcessed/furtherAnalysis/MACS2/flagstat.txt

        # echo $npeaks | awk -v ip=$ip -v OFS='\t' '{ print ip, $1 }' | cat $peak_count_header - > ${ip}_${peaktype}s.count_mqc.tsv

        # Maybe im too flexible, i could add -f 0.10 or -F 0.10 to ask at least for 10% overlap
        # awk -v FRiP=$FRiP -v OFS='\t' '{print ip, FRiP}' | cat $frip_score_header - > ${ip}_${peaktype}s.FRiP_mqc.tsv
    fi
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
                    ${extraAnnot} \
                    -cpu ${SLURM_CPUS_PER_TASK} \
                    -annStats ${outpath}/HOMER/peakAnnotation/${label}_${peaktype}.annotateStats.txt \
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
        peakfiles=$(find ${peaksPath}/${cell}_*peaks.${peaktype} | \
                    { grep -v -e "_input" -v -e "_IgG" || :; } | tr '\n' ',' )
        labels=$(find ${peaksPath}/${cell}_*peaks.${peaktype} -exec basename {} \; | \
                    { grep -v -e "_input" -v -e "_IgG" || :; })
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
                    { grep -v -e "_input" -v -e "_IgG" || :; } | tr '\n' ',' )
        labels=$(find ${annotPath}/${cell}_*${peaktype}.annotatePeaks.txt -exec basename {} \; | \
                    { grep -v -e "_input" -v -e "_IgG" || :; } )
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


chip="ATAC"
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
                { grep -v -e "_input" -v -e "_IgG" || :; } )
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




############################################################
# Annotate consensus peaks with HOMER, and add annotation to boolean output file
############################################################

if [ ! -e ${outpath}/HOMER/consensusPeaks ]; then
	mkdir -p ${outpath}/HOMER/consensusPeaks
fi

echo -e "Starting consensus peak annotations ------------------------------\n"


chip="ATAC"
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
                ${extraAnnot} \
                -cpu ${SLURM_CPUS_PER_TASK} \
                -annStats ${outpath}/HOMER/consensusPeaks/${prefix}.annotateStats.txt \
                > ${outpath}/HOMER/consensusPeaks/${prefix}.annotatePeaks.txt

        cut -f2- ${outpath}/HOMER/consensusPeaks/${prefix}.annotatePeaks.txt | \
            awk 'NR==1; NR > 1 {print $0 | "sort -T '.' -k1,1 -k2,2n"}' | \
            cut -f6- > ${outpath}/HOMER/consensusPeaks/tmp.txt
        paste ${outpath}/peakCalling/MACS2/consensusPeaks/${prefix}.boolean.txt \
            ${outpath}/HOMER/consensusPeaks/tmp.txt > ${boolAnotMatr}

        # We add a column for annotations regarding repeat elements
        if [[ $repeatsPath != "FALSE" ]]; then
            annotatePeaks.pl \
                    ${consensusPeakBed} \
                    ${speciesGenome} \
                    -gid \
                    -ann ${repeatsPath} \
                    -cpu ${SLURM_CPUS_PER_TASK} \
                    > ${outpath}/HOMER/consensusPeaks/${prefix}.annotatePeaks_rep.txt

            cut -f2- ${outpath}/HOMER/consensusPeaks/${prefix}.annotatePeaks_rep.txt | \
                awk 'NR==1; NR > 1 {print $0 | "sort -T '.' -k1,1 -k2,2n"}' | \
                cut -f7 > ${outpath}/HOMER/consensusPeaks/tmp.txt
            rm ${outpath}/HOMER/consensusPeaks/${prefix}.annotatePeaks_rep.txt
            # rename header and paste to consensus table
            sed -i 's/Annotation/Repeats Annotation/g' ${outpath}/HOMER/consensusPeaks/tmp.txt
            paste ${boolAnotMatr} \
                ${outpath}/HOMER/consensusPeaks/tmp.txt > ${outpath}/HOMER/consensusPeaks/tmp2.txt
            mv ${outpath}/HOMER/consensusPeaks/tmp2.txt ${boolAnotMatr}
        fi
    fi
done

echo -e "consensus peak annotations - Finished ---------------------\n"


###########################################################
# Count reads in consensus peaks with featureCounts
##########################################################

featureCpath=${outpath}/peakCalling/MACS2/consensusPeaks/featureCounts
if [ ! -e ${featureCpath} ]; then
	mkdir -p ${featureCpath}
fi
cd ${featureCpath}

echo -e "Starting consensus featureCounts -----------------------\n"

chip="ATAC"
for peaktype in narrowPeak broadPeak; do
    prefix="${chip}_${peaktype}_consensusPeaks"
    peakFiles=$(find -L ${outpath}/peakCalling/MACS2/peaks/*.${peaktype} -maxdepth 1  -type f ! -size 0 | \
                { grep -v -e "_input" -v -e "_IgG" || :; } )
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


###########################################################
# Read count to CPM
###########################################################

echo -e "Starting consensus CPM -----------------------\n"
chip="ATAC"
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
        sname=`basename $bam  | sed 's/.sort.rmdup.*.bam//g'`
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


###########################################################
# Make final table with all data
##########################################################

echo -e "Starting final merge table -----------------------\n"

# takes almost no time, so will repeat it always
python ${scriptsPath}/ChIP/cluster/02_NR_gatherAllInTable.py \
            ${outpath}/HOMER/consensusPeaks \
            ${outpath}/peakCalling/MACS2/consensusPeaks/featureCounts \
            ${outpath}/peakCalling/MACS2/consensusPeaks/annotWholeTable

echo -e "Final merge table- Finished ----------------------\n"



exit 0
