#!/bin/bash
## Initial SBATCH commands
#SBATCH --job-name=bcl2fastq
#SBATCH --time=07:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=12
#SBATCH -p short
#SBATCH -o /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.out
#SBATCH -e /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.err

## VARIABLES TO CHANGE
outdir=/home/jmendietaes/data/2021/singleCell/sequencedData/210806_VH00461_59_AAALHW2M5/
indir=/datos/intercambio/eguruce/210806_VH00461_59_AAALHW2M5_2
nSample=4
sampleSheet=/home/jmendietaes/data/2021/singleCell/sequencedData/demultiplex/Datasheet_CRISPRScreening.csv



## Load Software via modules/spack

module load bcl2fastq2/2.20.0-foss-2018b

## Run script
# Number of threads to write FASTQ data must be lower than number of samples
if [[ ${SLURM_CPUS_PER_TASK} -gt ${nSample} ]]; then
  wcpu=${nSample}
else
  wcpu=${SLURM_CPUS_PER_TASK}
fi


mkdir -p $outdir

bcl2fastq --runfolder-dir ${indir} \
        --output-dir  ${outdir}  \
        --sample-sheet $sampleSheet \
        --use-bases-mask Y28,I10,I10,Y84 \
        --minimum-trimmed-read-length 50 \
        --mask-short-adapter-read 7 \
        --no-lane-splitting \
        --barcode-mismatches 1 \
        -r ${SLURM_CPUS_PER_TASK} \
        -p ${SLURM_CPUS_PER_TASK} \
        -w ${wcpu} \
        -l INFO >> Ainhoa.dmux.MM0.log.txt 2>&1

