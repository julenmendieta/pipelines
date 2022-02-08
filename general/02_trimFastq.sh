#!/bin/bash
# -*- ENCODING: UTF-8 -*-

##===============================================================================
## SLURM VARIABLES
#SBATCH --job-name=trim
#SBATCH --cpus-per-task=8
#SBATCH --mem=5Gb
#SBATCH --time=1-00:00:00
#SBATCH -p medium
#SBATCH -o /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.out  
#SBATCH -e /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.err 

# HOW TO RUN ME
#sbatch /home/jmendietaes/programas/PhD/general/02_trimFastq.sh

# length for read1 and read2
r1Length=28
r2Length=84
# paths
mainpath="/home/jmendietaes/data/2021/singleCell/sequencedData/dataByExperiment"
insidemp="inVivo_OP2_ckit_14d_1"



# Number of threads
N=${SLURM_CPUS_PER_TASK}

# modules
module load Trimmomatic/0.38-Java-1.8
#module load cutadapt/1.18-foss-2018b-Python-2.7.15

for m in $insidemp; do
    cd ${mainpath}/${m}
    for novo in ${mainpath}/${m}/Novo*/gRNA; do
        read1=$(find ${novo}/*_R1_*gz)
        read1N=$(basename $read1)
        read1Dir=$(dirname $read1)
        read2=$(find ${novo}/*_R2_*gz)
        read2N=$(basename $read2)
        read2Dir=$(dirname $read2)

        #cutadapt -l ${r1Length} -o "${read1Dir}/trim_${read1N}" ${read1}
        java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.38.jar SE -threads $N ${read1} "${read1Dir}/trim_${read1N}" CROP:${r1Length}
        #cutadapt -l ${r2Length} -o "${read2Dir}/trim_${read2N}" ${read2}
        java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.38.jar SE -threads $N ${read2} "${read2Dir}/trim_${read2N}" CROP:${r2Length}

    done

    for novo in ${mainpath}/${m}/Novo*/mRNA; do
        read1=$(find ${novo}/*_R1_*gz)
        read1N=$(basename $read1)
        read1Dir=$(dirname $read1)
        read2=$(find ${novo}/*_R2_*gz)
        read2N=$(basename $read2)
        read2Dir=$(dirname $read2)

        #cutadapt -l ${r1Length} -o "${read1Dir}/trim_${read1N}" ${read1}
        #cutadapt -l ${r2Length} -o "${read2Dir}/trim_${read2N}" ${read2}
        java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.38.jar SE -threads $N ${read1} "${read1Dir}/trim_${read1N}" CROP:${r1Length}
        java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.38.jar SE -threads $N ${read2} "${read2Dir}/trim_${read2N}" CROP:${r2Length}


    done
done
