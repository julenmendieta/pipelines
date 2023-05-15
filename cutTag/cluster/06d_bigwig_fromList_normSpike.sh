#!/bin/bash
# -*- ENCODING: UTF-8 -*-

##===============================================================================
## SLURM VARIABLES
#SBATCH --job-name=specificBigwig
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH -p short
#SBATCH -o /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.out  
#SBATCH -e /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.err 

# HOW TO RUN ME
# add in the list bellow the bam file names (but the .bam termination)
#sbatch /home/jmendietaes/programas/pipelines/cutTag/cluster/06d_bigwig_fromList_normSpike.sh

#files="LSK-Chd4_ATAC5-merged.sort.rmdup.rmblackls.rmchr.Tn5 \
#LSK-Men1_ATAC5-merged.sort.rmdup.rmblackls.rmchr.Tn5"

# Path to base name to point to .fa and .EpiListRenamed
# .fa would be the fasta file with the spike-in sequences
# .EpiListRenamed contains a list of in-house ids for epigenetic marks
#   that are included at CUTANA_K-MetStat
spikeInRef="/home/jmendietaes/referenceGenomes/spikeIn/SNAP-CUTANA_K-MetStat"



files="DM-Kmt2a_H3K4me1_CUTTAG2_S9.sort.rmUnM.q30.rmchr.Tn5 \
DM-Kmt2a_H3K4me3_CUTTAG2_S10.sort.rmUnM.q30.rmchr.Tn5 \
DM-Kmt2a_K27ac_CUTTAG3_S12.sort.rmUnM.q30.rmchr.Tn5 \
DM-Kmt2b-As_H3K4me1_CUTTAG4_S12.sort.rmUnM.q30.rmchr.Tn5 \
DM-Kmt2b-As_H3K4me3_CUTTAG4_S13.sort.rmUnM.q30.rmchr.Tn5 \
DM-Kmt2b-Br_H3K4me1_CUTTAG4_S14.sort.rmUnM.q30.rmchr.Tn5 \
DM-Kmt2b-Br_H3K4me3_CUTTAG4_S15.sort.rmUnM.q30.rmchr.Tn5 \
DM-Kmt2d_H3K4me1_CUTTAG2_S13.sort.rmUnM.q30.rmchr.Tn5 \
DM-Kmt2d_H3K4me3_CUTTAG2_S7.sort.rmUnM.q30.rmchr.Tn5 \
DM-Kmt2d_IgG_CUTTAG2_S8.sort.rmdup.q30.rmchr.Tn5 \
DM-Kmt2d_K27ac_CUTTAG3_S13.sort.rmUnM.q30.rmchr.Tn5 \
DM-NTC_H3K4me1_CUTTAG2_S14.sort.rmUnM.q30.rmchr.Tn5 \
DM-NTC_H3K4me1_CUTTAG4_S20.sort.rmUnM.q30.rmchr.Tn5 \
DM-NTC_H3K4me3_CUTTAG2_S11.sort.rmUnM.q30.rmchr.Tn5 \
DM-NTC_H3K4me3_CUTTAG4_S21.sort.rmUnM.q30.rmchr.Tn5 \
DM-NTC_IgG_CUTTAG2_S12.sort.rmdup.q30.rmchr.Tn5 \
DM-NTC_K27ac_CUTTAG3_S14.sort.rmUnM.q30.rmchr.Tn5 \
DM-Setd1a_H3K4me1_CUTTAG4_S16.sort.rmUnM.q30.rmchr.Tn5 \
DM-Setd1a_H3K4me3_CUTTAG4_S17.sort.rmUnM.q30.rmchr.Tn5 \
DM-Stat5a_K27ac_CUTTAG3_S15.sort.rmUnM.q30.rmchr.Tn5 \
DM-Wdr82_H3K4me1_CUTTAG4_S18.sort.rmUnM.q30.rmchr.Tn5 \
DM-Wdr82_H3K4me3_CUTTAG4_S19.sort.rmUnM.q30.rmchr.Tn5 \" # \
#GMP-Brd9_ATACTraspl-merged.sort.rmdup.rmblackls.rmchr.Tn5"

basePath="/home/jmendietaes/data/2021/cutTag/allProcessed"
bamBase="${basePath}/bamfiles/valid/03_properAnalysis"


REFERENCE_DIR="/home/jmendietaes/referenceGenomes/mm10_reordered/mm10.reordered"
chr_genome_size=$REFERENCE_DIR".sizes"
wigToBigWig="/home/jmendietaes/programas/binPath/wigToBigWig"


echo -e "Starting BigWigs normalization---------------------------\n"

# check content of eleventh line of step control file
for fi in ${files}; do
    echo ${fi}

    bamPath="${bamBase}/${fi}.bam"
    bigWigOut2="${basePath}/BigWig/valid/${fi}.spike.bw"

    # Check if we can use the spike-in norm
    filename=$( echo $fi | sed 's/.sort.rmUnM.*//g')
    epiName=(${filename//_/ }); epiName=${epiName[1]}; 
    epiName=(${epiName//-/ }); epiName=${epiName[0]}
    summaryFile="${basePath}/QC/summary_${filename}.txt"

    if grep -v "#" ${spikeInRef}.EpiListRenamed | grep -Fxq "${epiName}"; then
        # All the information is at ${summaryFile}
        # get read1 counts and multiply by 2
        allRead=$(grep -A 2 "READ COUNTS" ${summaryFile} | tail -n 1)
        allRead=(${allRead// / }); allRead=${allRead[2]}; 
        allRead=$(($allRead*2))

        # get % of reads aligning to reference
        genomPerce=$(grep -A 7 "SAMTOOLS FLAGSTAT - DUPLICATES" ${summaryFile} | \
                    tail -n 1)
        genomPerce=(${genomPerce//(/ }); genomPerce=${genomPerce[4]}; 
        genomPerce=$(echo ${genomPerce} | sed 's/%//g')

        ## get % of reads aligning to spike-in
        # Get original chip name and look for counts to its spike-in
        # column 1 original name, column 2 in-house chip id
        epiOrigName=$(grep ${epiName} ${spikeInRef}.EpiListConversion | \
                        awk '{print $1}')
        lineNs=$(grep -n ${epiOrigName} ${spikeInRef}.fa | \
                cut --delimiter=":" --fields=1)
        # remove headers from count
        lineNs=$(for li in ${lineNs}; do echo $(((${li}+1)/2)); done)
        spikeR1=$(for li in ${lineNs}; do 
                sed "${li}q;d" ${basePath}/spikeQC/counts/${filename}_R1.txt;
            done | awk '{sum+=$1;} END{print sum;}')
        spikeR2=$(for li in ${lineNs}; do 
                sed "${li}q;d" ${basePath}/spikeQC/counts/${filename}_R2.txt;
            done | awk '{sum+=$1;} END{print sum;}')
        spikeAdd=$((${spikeR1} + ${spikeR2}))
        spikePerce=$(echo "print((${spikeAdd}/${allRead})*100)"  | python)
        scale_factor=$(echo "print(${genomPerce}/${spikePerce})"  | python)

        bamCoverage --binSize 1 --scaleFactor ${scale_factor} \
            -b ${bamPath} -of bigwig \
            -o ${bigWigOut2} --numberOfProcessors $SLURM_CPUS_PER_TASK



    else
        echo "ERROR: File cannot be scaled with this set of spike-ins"
        echo $fi
        exit 1
    fi


    echo -e "BigWigs second - done ---------------------------------------------\n"
done  
#sortBed -i "$i" | genomeCoverageBed -bg -i stdin -g ~/alignments/chr_genome_size/chr_hg38_size/hg38.chrom.sizes | wigToBigWig -clip stdin ~/alignments/chr_genome_size/chr_hg38_size/hg38.chrom.sizes ${f%.*PE*}.bw


echo -e "END --------------------------------------------------"
