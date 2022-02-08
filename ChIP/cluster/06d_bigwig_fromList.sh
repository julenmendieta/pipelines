#!/bin/bash
# -*- ENCODING: UTF-8 -*-

##===============================================================================
## SLURM VARIABLES
#SBATCH --job-name=peakAnalysis
#SBATCH --cpus-per-task=12
#SBATCH --mem=12G
#SBATCH --time=08:00:00
#SBATCH -p short
#SBATCH -o /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.out  
#SBATCH -e /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.err 

# HOW TO RUN ME
# add in the list bellow the bam file names (but the .bam termination)
#sbatch /home/jmendietaes/programas/PhD/ChIP/cluster/06d_bigwig_fromList.sh

files="MEP_Brd9-merged-sub105350855.sort.rmdup.rmblackls.rmchr \
MEP_Kmt2a-sub56264160_ChIP11_S7.sort.rmdup.rmblackls.rmchr \
MEP_Kmt2d-sub43061841_ChIP11_S3.sort.rmdup.rmblackls.rmchr \
GMP_Brd9-merged-sub105350855.sort.rmdup.rmblackls.rmchr \
GMP_Kmt2a-sub56264160_ChIP11_S9.sort.rmdup.rmblackls.rmchr \
GMP_Kmt2d-sub43061841_ChIP11_S4.sort.rmdup.rmblackls.rmchr \
GMP_Smarcb1_ChIP11_S1.sort.rmdup.rmblackls.rmchr"

basePath="/home/jmendietaes/data/2021/chip/allProcessed"
REFERENCE_DIR="/home/jmendietaes/referenceGenomes/mm10_reordered/mm10.reordered"
chr_genome_size=$REFERENCE_DIR".sizes"
wigToBigWig="/home/jmendietaes/programas/binPath/wigToBigWig"


echo -e "Starting BigWigs normalization---------------------------\n"

# check content of eleventh line of step control file
for fi in ${files}; do
    echo ${fi}
    bamPath="${basePath}/bamfiles/valid/subsampled_noIgG/${fi}.bam"
    bigWigOut2="${basePath}/BigWig/valid/${fi}.norm.bw"
    bamCoverage --binSize 5 --normalizeUsing CPM --exactScaling \
    -b ${bamPath} -of bigwig \
    -o ${bigWigOut2} --numberOfProcessors $SLURM_CPUS_PER_TASK

    echo -e "BigWigs second - done ---------------------------------------------\n"
done  
#sortBed -i "$i" | genomeCoverageBed -bg -i stdin -g ~/alignments/chr_genome_size/chr_hg38_size/hg38.chrom.sizes | wigToBigWig -clip stdin ~/alignments/chr_genome_size/chr_hg38_size/hg38.chrom.sizes ${f%.*PE*}.bw


echo -e "END --------------------------------------------------"