#!/bin/bash
# -*- ENCODING: UTF-8 -*-

##===============================================================================
## SLURM VARIABLES
#SBATCH --job-name=specificBigwig
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH -p short
#SBATCH -o /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.out  
#SBATCH -e /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.err 

# HOW TO RUN ME
# add in the list bellow the bam file names (but the .bam termination)
#sbatch /home/jmendietaes/programas/pipelines/mRNA/cluster/06d2_bigwig_fromList_TE.sh

#files="LSK-Chd4_ATAC5-merged.sort.rmdup.rmblackls.rmchr.Tn5 \
#LSK-Men1_ATAC5-merged.sort.rmdup.rmblackls.rmchr.Tn5"

# Set to "yes" if you want to generate a Bigwig file with info from both strands
bothStrands="no"
# Set to "yes" if you want to generate a Bigwig file per strand
perStrand="yes"

files="DMd4-Trim28-1-merged.sort.rmdup.rmblackls.rmchr \
DMd7-NTC-1-merged.sort.rmdup.rmblackls.rmchr \
DMd7-Pu1-1-merged.sort.rmdup.rmblackls.rmchr \
DMd7-Trim28-1-merged.sort.rmdup.rmblackls.rmchr" # \
#GMP-Brd9_ATACTraspl-merged.sort.rmdup.rmblackls.rmchr.Tn5"

basePath="/home/jmendietaes/data/2021/mRNA/allProcessed/TEtranscripts"
bamBase="${basePath}/bamfiles/wholeGenome/valid/merged/03_TEtranscripts"

export PATH="/home/jmendietaes/programas/miniconda3/bin:$PATH"
REFERENCE_DIR="/home/jmendietaes/referenceGenomes/mm10_reordered/mm10.reordered"
chr_genome_size=$REFERENCE_DIR".sizes"
wigToBigWig="/home/jmendietaes/programas/binPath/wigToBigWig"
bamCoverageP="/home/jmendietaes/programas/miniconda3/bin/bamCoverage"

echo -e "Starting BigWigs normalization---------------------------\n"

# Get Bigwig with both strands
if [[ ${bothStrands} == "yes" ]]; then 
    mkdir -p ${basePath}/BigWig/stranded
    for fi in ${files}; do
        echo ${fi}
        bamPath="${bamBase}/${fi}.bam"
        bigWigOut2="${basePath}/BigWig/${fi}.norm.bw"
        
        ${bamCoverageP} --binSize 5 --normalizeUsing CPM --exactScaling \
        -b ${bamPath} -of bigwig \
        -o ${bigWigOut2} --numberOfProcessors $SLURM_CPUS_PER_TASK

        echo -e "BigWigs second - done ---------------------------------------------\n"
    done  
fi
#sortBed -i "$i" | genomeCoverageBed -bg -i stdin -g ~/alignments/chr_genome_size/chr_hg38_size/hg38.chrom.sizes | wigToBigWig -clip stdin ~/alignments/chr_genome_size/chr_hg38_size/hg38.chrom.sizes ${f%.*PE*}.bw

# Get Bigwig per strand
# https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html#separate-tracks-for-each-strand
# We use Illumina dUTP
# Rd1 maps antisense strand
# Rd2 maps sense strand
# The sense strand is the strand of DNA that has the same sequence as the mRNA, 
# which takes the antisense strand as its template during transcription
# Default way if where Read 2 (R2) is in the direction of RNA strand (reverse-stranded library)
if [[ ${perStrand} == "yes" ]]; then 
    for fi in ${files}; do
        echo ${fi}
        bamPath="${bamBase}/${fi}.bam"
        
        bigWigOut2="${basePath}/BigWig/stranded/${fi}.norm.forward.bw"
        ${bamCoverageP} --binSize 5 --normalizeUsing CPM --exactScaling \
        -b ${bamPath} -of bigwig \
        -o ${bigWigOut2} --numberOfProcessors $SLURM_CPUS_PER_TASK \
        --filterRNAstrand forward

        bigWigOut2="${basePath}/BigWig/stranded/${fi}.norm.reverse.bw"
        ${bamCoverageP} --binSize 5 --normalizeUsing CPM --exactScaling \
        -b ${bamPath} -of bigwig \
        -o ${bigWigOut2} --numberOfProcessors $SLURM_CPUS_PER_TASK \
        --filterRNAstrand reverse

        echo -e "BigWigs second - done ---------------------------------------------\n"
    done  
fi


echo -e "END --------------------------------------------------"
