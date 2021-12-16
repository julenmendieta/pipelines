#!/bin/bash
# -*- ENCODING: UTF-8 -*-

##===============================================================================
## SLURM VARIABLES
#SBATCH --job-name=hichip
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=10:00:00
#SBATCH -p short
#SBATCH -o /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.out  
#SBATCH -e /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.err 


# HOW TO RUN ME
#for i in *fastq.gz; do echo $i | sed 's/_R._001.fastq.gz//g' ; done | sort | uniq > samplesNames.txt  
# N=`cat samplesNames.txt | wc -l`
# sbatch --array=1-${N} /home/jmendietaes/programas/PhD/HiChIP/cluster/01_Dovetail_generateBAM.sh \
#/home/jmendietaes/data/2021/HiChIP/sequencedData/NextSeq2000.RUN25.20210923 \
#/home/jmendietaes/data/2021/HiChIP/allProcessed \
#/home/jmendietaes/referenceGenomes/mm10_reordered/mm10.reordered

# OBJECTIVE
# process hichip data as stated in https://hichip.readthedocs.io/en/latest


## GLOBAL VARIABLES
PROJECT_DIR=$1
#PROJECT_DIR="/home/jmendietaes/data/2021/HiChIP/sequencedData/NextSeq2000.RUN25.20210923"
# In my case, the fastq files are in project folder inside demux_fastq

RAW_FASTQ_DIR=$PROJECT_DIR"/demux_fastq"
# not used, but i will store here final files
#RESULTS_DIR=$PROJECT_DIR"/pipelineOut/results"
EDITED_DIR=$PROJECT_DIR"/pipelineOut"
FASTQ_DIR=$EDITED_DIR"/fastq" 
#libMetricsP="${EDITED_DIR}/libMetrics"

## path to the main fully processed data output
basePath=$2
#basePath="/home/jmendietaes/data/2021/HiChIP/allProcessed"
# here we will stored the final filtered bam files
bamsPath="${basePath}/bamfiles"


##  REFERENCE Genome info
REFERENCE_DIR=$3
#REFERENCE_DIR="/home/jmendietaes/referenceGenomes/mm10_reordered/mm10.reordered"
# REFERENCE_DIR includes the path and naming base of the reference genome, we will
#   add .sizes, .blacklist.bed to the base to get the rest of paths. Genome
#   indexes must follow the same base name

GenomeIndex=$REFERENCE_DIR
#chrOrder=$REFERENCE_DIR"/mm10_Bowtie2/names.txt"
#chr_genome_size=$REFERENCE_DIR".sizes"
#BlackList=$REFERENCE_DIR".blacklist.bed"
#wigToBigWig="/home/jmendietaes/programas/binPath/wigToBigWig"
#picardPath='/home/jmendietaes/programas/picard/picard.jar'


## load modules
hichipScripts="/home/jmendietaes/programas/HiChiP"
export PATH="/home/jmendietaes/programas/miniconda3/bin:$PATH"
module load BWA/0.7.17-foss-2018b
module load SAMtools/1.12-GCC-10.2.0

##===============================================================================

# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

##===============================================================================

cd $EDITED_DIR

#Choose file
# If we always get in list all files from folder
#FILES=($(ls $RAW_FASTQ_DIR | grep "fastq.gz" | sed -r 's/.{16}$//' | uniq))
# If we want to do only the files in samplesNames.txt file. Perfecto for re-running 
#dead jobs
FILES=($(cat $RAW_FASTQ_DIR/samplesNames.txt))
filename=${FILES[$SLURM_ARRAY_TASK_ID - 1]}
echo -e $SLURM_JOB_NODELIST 
#https://slurm.schedmd.com/faq.html#task_prolog
echo "print =========================================="
echo "print SLURM_JOB_ID = $SLURM_JOB_ID"
echo "print SLURM_JOB_NODELIST = $SLURM_JOB_NODELIST"
echo "print =========================================="
echo $filename
##===============================================================================


# get some paths
read1_path="${RAW_FASTQ_DIR}/${filename}_R1_001.fastq.gz"
read2_path="${RAW_FASTQ_DIR}/${filename}_R2_001.fastq.gz"
tempdir="${EDITED_DIR}/tmp"
statsOut="${EDITED_DIR}/QC/pairtools_${filename}"
outpair="${EDITED_DIR}/pairtools/outpairs.gz"
outbam="${bamsPath}/${filename}.PT.bam"
samFile="${EDITED_DIR}/BAM/${filename}.sam"


if [ ! -e ${EDITED_DIR} ]; then
    mkdir -p ${EDITED_DIR}
    mkdir -p ${EDITED_DIR}/QC
    mkdir -p ${EDITED_DIR}/pairtools
    mkdir -p ${tempdir}
    mkdir -p ${EDITED_DIR}/BAM
fi


SLURM_CPUS_PER_TASK=16

## fast version
# here we will:
#- Align
#- Record valid ligation events
#- sort pairsam file
#- remove PCR duplicates
#- generate pairs and bam files
bwa mem -5SP -T0 -t${SLURM_CPUS_PER_TASK} ${REFERENCE_DIR}.fa ${read1_path} ${read2_path}| \
pairtools parse --min-mapq 40 --walks-policy 5unique \
--max-inter-align-gap 30 --nproc-in ${SLURM_CPUS_PER_TASK} --nproc-out ${SLURM_CPUS_PER_TASK} --chroms-path ${REFERENCE_DIR}.fa | \
pairtools sort --tmpdir=${tempdir} --nproc ${SLURM_CPUS_PER_TASK} |pairtools dedup --nproc-in ${SLURM_CPUS_PER_TASK} \
--nproc-out ${SLURM_CPUS_PER_TASK} --mark-dups --output-stats "${statsOut}.txt" | pairtools split --nproc-in ${SLURM_CPUS_PER_TASK} \
--nproc-out ${SLURM_CPUS_PER_TASK} --output-pairs ${outpair} --output-sam - | samtools view -bS -@${SLURM_CPUS_PER_TASK} | \
samtools sort -@${SLURM_CPUS_PER_TASK} -o ${outbam} ;samtools index ${outbam}

## slow version
# # Alignment
# bwa mem -5SP -T0 -t${SLURM_CPUS_PER_TASK} ${REFERENCE_DIR}.fa \
#                 ${read1_path} ${read2_path} -o ${samFile}

# # Recording valid ligation events
# parsedSam="${EDITED_DIR}/BAM/${filename}.parsed.sam"
# pairtools parse --min-mapq 40 --walks-policy 5unique \
#         --max-inter-align-gap 30 --nproc-in ${SLURM_CPUS_PER_TASK} \
#         --nproc-out ${SLURM_CPUS_PER_TASK} \
#         --chroms-path ${REFERENCE_DIR}.fa ${samFile} > ${parsedSam}

# # sorting pairsam file
# sortedSam="${EDITED_DIR}/BAM/${filename}.parsed.sorted.sam"
# pairtools sort --nproc ${SLURM_CPUS_PER_TASK} --tmpdir=${tempdir} ${parsedSam} > ${sortedSam}

# # Remove PCR duplicates
# dedupSAM="${EDITED_DIR}/BAM/${filename}.parsed.sorted.dedup.sam"
# pairtools dedup --nproc-in ${SLURM_CPUS_PER_TASK} \
#             --nproc-out ${SLURM_CPUS_PER_TASK} --mark-dups \
#             --output-stats "${statsOut}.txt" --output ${dedupSAM} \
#              ${sortedSam}

# # Generate pairs and bam files
# unsortSAM="${EDITED_DIR}/BAM/${filename}.parsed.sorted.dedup.unsort.sam"
# pairtools split --nproc-in ${SLURM_CPUS_PER_TASK} \
#         --nproc-out ${SLURM_CPUS_PER_TASK} --output-pairs ${outpair} \
#          --output-sam ${unsortSAM} ${dedupSAM}

# samtools sort -@${SLURM_CPUS_PER_TASK} -T ${tempdir}/tempfile.bam -o ${outbam} ${unsortSAM}
# samtools index ${outbam}

python3 ${hichipScripts}/get_qc.py -p "${statsOut}.txt" > "${statsOut}_processed.txt"

peaksBed='/home/jmendietaes/data/2021/chip/allProcessed/furtherAnalysis/subsampled_noIgG/peakCalling/MACS2/peaks/Mye_Smarcb1-merged-sub173331792_peaks.narrowPeak'

${hichipScripts}/enrichment_stats.sh -g ${REFERENCE_DIR}.fa -b ${outbam} -p ${peaksBed} -t ${SLURM_CPUS_PER_TASK} -x ${statsOut}


python3 ${hichipScripts}/plot_chip_enrichment.py -bam ${outbam} -peaks ${peaksBed} -output ${statsOut}_enrichment.png>
