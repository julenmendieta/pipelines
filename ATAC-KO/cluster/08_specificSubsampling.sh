#!/bin/bash
# -*- ENCODING: UTF-8 -*-

##===============================================================================
## SLURM VARIABLES
#SBATCH --job-name=specificSumbsampling
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH -p short
#SBATCH -o /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.out  
#SBATCH -e /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.err 


# HOW TO RUN ME
#sbatch /home/jmendietaes/programas/pipelines/ATAC-KO/cluster/08_specificSubsampling.sh

# OBJECTIVE
# call peaks, annotate, make consensus peaks, get reads in peaks, CPM and merge
# all in same table

module load Sambamba/0.7.0
# function to join elements from array
function join_by { local IFS="$1"; shift; echo "$*"; }


## INPUT PARAMETERS
# bam input and output paths
inpath=/home/jmendietaes/data/2021/ATAC/allProcessed/bamfiles/valid/mergedReplicates/08e_koATAC
outpath=/home/jmendietaes/data/2021/ATAC/allProcessed/bamfiles/valid/mergedReplicates/08e_koATAC_subsampled
# Cell name and batch to subsample
batch=ATAC7
cell=DM
# Counted as R1 + R2
minRead=28369833
# Number of CPU to use
#SLURM_CPUS_PER_TASK=8


################################################################################
# Version for merged files
#allbams=$(find ${inpath}/${cell}-*${batch}-*bam -printf "%f\n" | \
#            { grep -v -e "-Kat" || :; })
# Version for replicates
allbams=$(find ${inpath}/${cell}-*${batch}_*bam -printf "%f\n" | \
            { grep -v -e "-Kat" || :; })

for bam in ${allbams}; do
    nReads=$(samtools idxstats ${inpath}/${bam} | cut -f 3 | \
            awk '{sum+=$1;} END{print sum;}')
    fractionOfReads=$(echo "print(${minRead}/${nReads})" | python3)

    # If to be subsampled
    if (( $(echo $fractionOfReads'<'1 |bc -l) )) ; then
        echo "$bam"
        echo "Subtracted"
        mapLib=(${bam//\./ }); 
        # We add as approx the number of reads we wanted
        mapLib[0]=${mapLib[0]}"-sub${minRead}"
        subsName=$(join_by . ${mapLib[@]})

        sambamba view -h -t $SLURM_CPUS_PER_TASK -s $fractionOfReads \
                            -f bam --subsampling-seed=12345 ${inpath}/${bam} \
                            -o ${outpath}/${subsName}

        # We want final number of reads
        nReads=$(samtools idxstats ${outpath}/${subsName} | cut -f 3 | \
            awk '{sum+=$1;} END{print sum;}')
        mapLib=(${bam//\./ }); 
        # We add as approx the number of reads we wanted
        mapLib[0]=${mapLib[0]}"-sub${nReads}"
        subsName2=$(join_by . ${mapLib[@]})
        mv ${outpath}/${subsName} ${outpath}/${subsName2}
        mv ${outpath}/${subsName}.bai ${outpath}/${subsName2}.bai

    else
        echo "$bam"
        mapLib=(${bam//\./ }); 
        mapLib[0]=${mapLib[0]}"-noSub${nReads}"
        subsName=$(join_by . ${mapLib[@]})
        echo "Not subtracted: fraction=${fractionOfReads}"
        ln -s  ${inpath}/${bam} ${outpath}/${subsName}
        ln -s  ${inpath}/${bam}.bai ${outpath}/${subsName}.bai
    fi
    echo "Done"

done

# Get NTCs
# allbams=$(find ${inpath}/${cell}-*${batch}-*bam -printf "%f\n" | \
#             { grep -e "NTC" || :; })
# for bam in ${allbams}; do
#     nReads=$(samtools idxstats ${inpath}/${bam} | cut -f 3 | \
#             awk '{sum+=$1;} END{print sum;}')
#     echo "$bam"
#     mapLib=(${bam//\./ }); 
#     mapLib[0]=${mapLib[0]}"-noSub${nReads}"
#     subsName=$(join_by . ${mapLib[@]})

#     ln -s  ${inpath}/${bam} ${outpath}/${subsName}
#     ln -s  ${inpath}/${bam}.bai ${outpath}/${subsName}.bai
# done

# # Check average reads in NTC merged
# allbams=$(find ${inpath}/${cell}-*${batch}-*bam -printf "%f\n" | \
#             { grep -e "NTC" || :; })
# for bam in ${allbams}; do
#     nReads=$(samtools idxstats ${inpath}/${bam} | cut -f 3 | \
#             awk '{sum+=$1;} END{print sum;}')
#     echo "$bam  ${nReads}"
# done