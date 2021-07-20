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



#############################
#PEAK CALLING: MACS2
#############################


echo -e "Starting Peak calling -----------------------------------------------\n"


# Create output dir
if [ ! -e ${outpath}/peakCalling/MACS2/logs ]; then
    mkdir -p ${outpath}/peakCalling/MACS2/logs
    mkdir -p ${outpath}/peakCalling/MACS2/peaks
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
    cellControls=$(echo ${allbams} | tr ' ' '\n' | grep "${cell}_" | grep -e input -e IgG)
    label=$(basename $bam | cut -d '.' -f 1 | cut -d '_' -f 1,2,3)
    # for now we use the first control to normalise signal (by sort should be IgG)
    controlbam=`echo ${cellControls} | cut -d ' ' -f 1`

    total_reads=$(samtools view -c ${bam})

    # narrow peaks
    peaktype='narrowPeak'

    if [ ! -e ${outpath}/peakCalling/MACS2/peaks/${label}_peaks_${peaktype}.xls ]; then

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

    if [ ! -e ${outpath}/peakCalling/MACS2/peaks/${label}_peaks_${peaktype}.xls ]; then
        macs2 callpeak \
                -t ${bam} \
                -c ${controlbam} \
                --broad \
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


###############################
# HOMER: annotate peaks
##############################

# anotate peaks HOMER
# http://homer.ucsd.edu/homer/ngs/annotation.html

# Create output dir
if [ ! -e ${outpath}/HOMER ]; then
    mkdir -p ${outpath}/HOMER
fi


for bam in allbams; do
    for peaktype in 'narrowPeak' 'broadPeak'; do
        label=$(basename $bam | cut -d '.' -f 1 | cut -d '_' -f 1,2,3; done)
        annotatePeaks.pl \
                ${outpath}/peakCalling/MACS2/peaks/${label}_peaks.${peaktype} \
                ${speciesGenome} \
                -gid \
                -cpu ${SLURM_CPUS_PER_TASK} \
                > ${outpath}/HOMER/${label}_${peaktype}.annotatePeaks.txt
    done
done

###############################
# QC: MACS2 quality check plots
##############################
# plot metrics
if [ ! -e ${outpath}/peakCalling/MACS2/QC ]; then
    mkdir -p ${outpath}/peakCalling/MACS2/QC
fi

allCell=`for filename in ${allLabels}; do 
            mapLib=(${filename//_/ }); 
            echo ${mapLib[0]}; done | sort | uniq`

for cell in ${allCell}; do
    for peaktype in 'narrowPeak' 'broadPeak'; do
    
        # get comma separated string with all the peak files to compare
        peakfiles=$(find ${outpath}/peakCalling/MACS2/peaks/${cell}_*peaks.${peaktype} -printf "${outpath}/peakCalling/MACS2/peaks/%f,")
        labels=$(find ${outpath}/peakCalling/MACS2/peaks/${cell}_*peaks.${peaktype} -printf "%f ")
        labels=$(for la in $labels; do echo $la | sed "s/_peaks.${peaktype}//g"; done | tr '\n' ',' )

        # https://github.com/nf-core/chipseq
        Rscript ${scriptsPath}/ChIP/cluster/02_NR_plot_macs_qc.R \
                -i ${peakfiles} \
                -s ${labels} \
                -o ${outpath}/peakCalling/MACS2/QC \
                -p macs_peak # out prefix


        # get comma separated string with all the peak annotation files to compare
        annotfiles=$(find ${outpath}/HOMER/${cell}_*${peaktype}.annotatePeaks.txt -printf "${outpath}/HOMER/%f,")
        labels=$(find ${outpath}/HOMER/${cell}_*${peaktype}.annotatePeaks.txt -printf "%f ")
        labels=$(for la in $labels; do echo $la | sed "s/_${peaktype}.annotatePeaks.txt//g"; done | tr '\n' ',' )

        Rscript ${scriptsPath}/ChIP/cluster/02_NR_plot_homer_annotatepeaks.R \
                -i ${annotfiles} \
                -s ${labels} \
                -o ${outpath}/peakCalling/MACS2/QC \
                -p macs_annotatePeaks
    done
done


#########################
# CONSENSUS PEAKS ANALYSIS: Consensus peaks across samples, create boolean 
#     filtering file, SAF file for featureCounts and UpSetR plot for intersection
#########################

## We first do it in the whole dataset as a trial, but in the future ill do it only between the biological replicates
if [ ! -e ${outpath}/PeakCalling/MACS2/consensusPeaks ]; then
	mkdir -p ${outpath}/PeakCalling/MACS2/consensusPeaks
fi

chip="allmerged"
for peaktype in narrowPeak broadPeak; do
    # hay que usar estos
    if ${peaktype} == "narrowPeak"; then
        mergecols=`seq 2 10 | tr '\n' ','`
        expandparam='--is_narrow_peak'
    elif ${peaktype} == "broadPeak"; then
        mergecols=`seq 2 9 | tr '\n' ','`
        expandparam=''
    fi

    peakFiles=$(ls ${outpath}/peakCalling/MACS2/peaks/*.${peaktype})
    fileLabels=$(for f in $peakFiles; do echo $f |sed "s/_peaks.${peaktype}//g"; done | tr '\n' ',')


    prefix="${chip}_${peaktype}_consensusPeaks"

    sort -T '.' -k1,1 -k2,2n ${peakFiles} \
        | mergeBed -c $mergecols -o collapse > ${outpath}/peakCalling/MACS2/consensusPeaks/${prefix}.txt

    python ${scriptsPath}/ChIP/cluster/02_NR_macs2_merged_expand.py ${outpath}/peakCalling/MACS2/consensusPeaks/${prefix}.txt \
        ${fileLabels} \
        ${outpath}/peakCalling/MACS2/consensusPeaks/${prefix}.boolean.txt \
        $expandparam

    awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print $1, $2, $3, $4, "0", "+" }' \
        ${outpath}/peakCalling/MACS2/consensusPeaks/${prefix}.boolean.txt > ${outpath}/peakCalling/MACS2/consensusPeaks/${prefix}.bed
    echo -e "GeneID\tChr\tStart\tEnd\tStrand" > ${outpath}/peakCalling/MACS2/consensusPeaks/${prefix}.saf
    awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print $4, $1, $2, $3,  "+" }' \
        ${outpath}/peakCalling/MACS2/consensusPeaks/${prefix}.boolean.txt >> ${outpath}/peakCalling/MACS2/consensusPeaks/${prefix}.saf
    Rscript ${scriptsPath}/ChIP/cluster/02_NR_plot_peak_intersect.r \
        -i ${outpath}/peakCalling/MACS2/consensusPeaks/${prefix}.boolean.intersect.txt \
        -o ${outpath}/peakCalling/MACS2/consensusPeaks/${prefix}.boolean.intersect.plot.pdf
    
done

# Now we do it within the biological replicates 
# create file with ChIP names in bams folder
# for filename in *bam; do mapLib=(${filename//_/ }); mapLib=${mapLib[1]}; mapLib=(${mapLib//-/ }); echo ${mapLib[0]}; done | grep -v input | grep -v IgG | sort | uniq > chipMarcs.txt
# markFile="/home/jmendietaes/data/2021/chip/allProcessed/bamfiles/valid/chipMarcs.txt"
# PROTS=($(cat $bamsPath/chipMarcs.txt))

# if [ ! -e ${outpath}/PeakCalling/MACS2/consensusPeaks ]; then
# 	mkdir -p ${outpath}/PeakCalling/MACS2/consensusPeaks
# fi

# for peaktype in narrowPeak broadPeak; do
#     # hay que usar estos
#     if ${peaktype} == "narrowPeak"; then
#         mergecols=`seq 2 10 | tr '\n' ','`
#         expandparam='--is_narrow_peak'
#     elif ${peaktype} == "broadPeak"; then
#         mergecols=`seq 2 9 | tr '\n' ','`
#         expandparam=''
#     fi

#     for Ig_prot in ${PROTS}; do
#         protBams=$(find ${bamsPath}/*_${Ig_prot}*bam -printf "%f ")
#         protNames=$(for p in $protBams; do basename $p | sed 's/\.sort\.rmdup.*//g'; done)

#         peakFiles=$(for p in $protNames; do echo ${outpath}/MACS2/${p}_peaks.${peaktype}; done)
#         fileLabels=$(for f in $peakFiles; do echo $f |sed "s/_peaks.${peaktype}//g"; done | tr '\n' ',')


#         prefix = "${Ig_prot}.consensus_peaks"

#         sort -T '.' -k1,1 -k2,2n ${peakFiles} \
#             | mergeBed -c $mergecols -o collapse > ${outpath}/peakCalling/MACS2/consensusPeaks/${prefix}.txt

#         python ${scriptsPath}/ChIP/cluster/02_NR_macs2_merged_expand.py ${outpath}/peakCalling/MACS2/consensusPeaks/${prefix}.txt \
#             ${fileLabels} \
#             ${outpath}/peakCalling/MACS2/consensusPeaks/${prefix}.boolean.txt \
#             $expandparam

#     awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print $1, $2, $3, $4, "0", "+" }' \
#         ${outpath}/peakCalling/MACS2/consensusPeaks/${prefix}.boolean.txt > ${outpath}/peakCalling/MACS2/consensusPeaks/${prefix}.bed
#     echo -e "GeneID\tChr\tStart\tEnd\tStrand" > ${outpath}/peakCalling/MACS2/consensusPeaks/${prefix}.saf
#     awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print $4, $1, $2, $3,  "+" }' \
#         ${outpath}/peakCalling/MACS2/consensusPeaks/${prefix}.boolean.txt >> ${outpath}/peakCalling/MACS2/consensusPeaks/${prefix}.saf
#     Rscript ${scriptsPath}/ChIP/cluster/02_NR_plot_peak_intersect.r \
#         -i ${outpath}/peakCalling/MACS2/consensusPeaks/${prefix}.boolean.intersect.txt \
#         -o ${outpath}/peakCalling/MACS2/consensusPeaks/${prefix}.boolean.intersect.plot.pdf
    
# done
# Consensus peaks => R



############################################################
# Annotate consensus peaks with HOMER, and add annotation to boolean output file
############################################################

if [ ! -e ${outpath}/HOMER/consensusPeaks ]; then
	mkdir -p ${outpath}/HOMER/consensusPeaks
fi

chip="allmerged"
for peaktype in narrowPeak broadPeak; do

    prefix="${chip}_${peaktype}_consensusPeaks"

    annotatePeaks.pl \
            ${outpath}/peakCalling/MACS2/consensusPeaks/${prefix}.bed \
            ${speciesGenome} \
            -gid \
            -cpu ${SLURM_CPUS_PER_TASK} \
            > ${outpath}/HOMER/consensusPeaks/${prefix}.annotatePeaks.txt

    cut -f2- ${outpath}/HOMER/consensusPeaks/${prefix}.annotatePeaks.txt | \
        awk 'NR==1; NR > 1 {print $0 | "sort -T '.' -k1,1 -k2,2n"}' | \
        cut -f6- > ${outpath}/HOMER/consensusPeaks/tmp.txt
    paste ${outpath}/peakCalling/MACS2/consensusPeaks/${prefix}.boolean.txt \
        ${outpath}/HOMER/consensusPeaks/tmp.txt > ${outpath}/HOMER/consensusPeaks/${prefix}.boolean.annotatePeaks.txt
done
###########################################################
# Count reads in consensus peaks with featureCounts
##########################################################

# this only for the consensus between replicates, here no sense
# featureCounts \
#         -F SAF \
#         -O \
#         --fracOverlap 0.2 \
#         -T $task.cpus \
#         $pe_params \
#         -a $saf \
#         -o ${prefix}.featureCounts.txt \
#         ${bam_files.join(' ')}


###########################################################
# Differential analysis with DESeq2
###########################################################

# Esto parte del ${prefix}.featureCounts.txt del paso anterior
# featurecounts_deseq2.r \\
#         --featurecount_file $counts \\
#         --bam_suffix '$bam_ext' \\
#         --outdir ./ \\
#         --outprefix $prefix \\
#         --outsuffix '' \\
#         --cores $task.cpus \\
#         $vst
#     sed 's/deseq2_pca/deseq2_pca_${task.index}/g' <$deseq2_pca_header >tmp.txt
#     sed -i -e 's/DESeq2 /${antibody} DESeq2 /g' tmp.txt
#     cat tmp.txt ${prefix}.pca.vals.txt > ${prefix}.pca.vals_mqc.tsv
#     sed 's/deseq2_clustering/deseq2_clustering_${task.index}/g' <$deseq2_clustering_header >tmp.txt
#     sed -i -e 's/DESeq2 /${antibody} DESeq2 /g' tmp.txt
#     cat tmp.txt ${prefix}.sample.dists.txt > ${prefix}.sample.dists_mqc.tsv


###########################################################
# IGV
###########################################################


