#!/bin/bash
# -*- ENCODING: UTF-8 -*-

##===============================================================================
## SLURM VARIABLES
#SBATCH --job-name=specificBigwig
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH -p short
#SBATCH -o /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.out  
#SBATCH -e /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.err 

# HOW TO RUN ME
# add in the list bellow the bam file names (but the .bam termination)
#sbatch /home/jmendietaes/programas/pipelines/ATAC/cluster/06d_bigwig_fromList.sh

#files="LSK-Chd4_ATAC5-merged.sort.rmdup.rmblackls.rmchr.Tn5 \
#LSK-Men1_ATAC5-merged.sort.rmdup.rmblackls.rmchr.Tn5"

files="DMd4-Pu1-1_ATAC11-merged.sort.rmdup.rmblackls.rmchr.Tn5 \
DMd7-Pu1-1_ATAC12-merged.sort.rmdup.rmblackls.rmchr.Tn5 \
DMd7-Pu1-2_ATAC12-merged.sort.rmdup.rmblackls.rmchr.Tn5 \
DMd7-NTC-2_ATAC12-merged.sort.rmdup.rmblackls.rmchr.Tn5" # \
#GMP-Brd9_ATACTraspl-merged.sort.rmdup.rmblackls.rmchr.Tn5"

basePath="/home/jmendietaes/data/2021/ATAC/allProcessed"
bamBase="${basePath}/bamfiles/valid/08e_koATAC"


REFERENCE_DIR="/home/jmendietaes/referenceGenomes/mm10_reordered/mm10.reordered"
chr_genome_size=$REFERENCE_DIR".sizes"
wigToBigWig="/home/jmendietaes/programas/binPath/wigToBigWig"


echo -e "Starting BigWigs normalization---------------------------\n"

# check content of eleventh line of step control file
for fi in ${files}; do
    echo ${fi}
    bamPath="${bamBase}/${fi}.bam"
    bigWigOut2="${basePath}/BigWig/valid/${fi}.norm.bw"
    bamCoverage --binSize 5 --normalizeUsing CPM --exactScaling \
    -b ${bamPath} -of bigwig \
    -o ${bigWigOut2} --numberOfProcessors $SLURM_CPUS_PER_TASK

    echo -e "BigWigs second - done ---------------------------------------------\n"
done  
#sortBed -i "$i" | genomeCoverageBed -bg -i stdin -g ~/alignments/chr_genome_size/chr_hg38_size/hg38.chrom.sizes | wigToBigWig -clip stdin ~/alignments/chr_genome_size/chr_hg38_size/hg38.chrom.sizes ${f%.*PE*}.bw


echo -e "END --------------------------------------------------"
