#!/bin/bash
# -*- ENCODING: UTF-8 -*-

##===============================================================================
## SLURM VARIABLES
#SBATCH --job-name=microC
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --time=20:00:00
#SBATCH -p short
#SBATCH -o /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.out  
#SBATCH -e /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.err 


# HOW TO RUN ME
#for i in *fastq.gz; do echo $i | sed 's/_R._001.fastq.gz//g' ; done | sort | uniq > samplesNames.txt  
# N=`cat samplesNames.txt | wc -l`
# sbatch --array=1-${N} /home/jmendietaes/programas/pipelines/microC/cluster/01_Dovetail_generateBAM.sh \
#/home/jmendietaes/data/2021/microC/sequencedData/NextSeq2000.RUN137.230120 \
#/home/jmendietaes/data/2021/microC/allProcessed \
#/home/jmendietaes/referenceGenomes/mm10_reordered/mm10.reordered

# OBJECTIVE
# process microC data as stated in https://micro-c.readthedocs.io/en/latest/fastq_to_bam.html


## GLOBAL VARIABLES
PROJECT_DIR=$1
#PROJECT_DIR="/home/jmendietaes/data/2021/microC/sequencedData/NextSeq2000.RUN25.20210923"
# In my case, the fastq files are in project folder inside demux_fastq

RAW_FASTQ_DIR=$PROJECT_DIR"/demux_fastq"
# not used, but i will store here final files
#RESULTS_DIR=$PROJECT_DIR"/pipelineOut/results"
EDITED_DIR=$PROJECT_DIR"/pipelineOut"
FASTQ_DIR=$EDITED_DIR"/fastq" 
#libMetricsP="${EDITED_DIR}/libMetrics"

## path to the main fully processed data output
basePath=$2
#basePath="/home/jmendietaes/data/2021/microC/allProcessed"
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

# Path to Preseq binary for WC
preseq=/home/jmendietaes/programas/preseq_v2.0/preseq-3.1.2/installation/bin/preseq

## load modules
microCScripts="/home/jmendietaes/programas/pipelines/microC/cluster"
# set to lowercase yes or not
removeTemp="no"
export PATH="/home/jmendietaes/programas/miniconda3/bin:$PATH"
module load BWA/0.7.17-foss-2018b
module load SAMtools/1.12-GCC-10.2.0
module load CMake/3.21.1-GCCcore-11.2.0
module load HTSlib/1.11-GCC-10.2.0
##===============================================================================

# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

##===============================================================================

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


##===============================================================================

# get some paths
read1_path="${RAW_FASTQ_DIR}/${filename}_R1_001.fastq.gz"
read2_path="${RAW_FASTQ_DIR}/${filename}_R2_001.fastq.gz"
tempdir="${EDITED_DIR}/tmp"
statsOut="${basePath}/QC/pairtools_${filename}"
outpair="${basePath}/pairtools/${filename}_outpairs.gz"
samFile="${EDITED_DIR}/BAM/${filename}.sam"


if [ ! -e ${EDITED_DIR} ]; then
    mkdir -p ${EDITED_DIR}
    mkdir -p ${basePath}/QC
    mkdir -p ${basePath}/pairtools
    mkdir -p ${tempdir}
    mkdir -p ${EDITED_DIR}/BAM
fi

cd $EDITED_DIR

#############
# Alignment and duplicates removal
#############

## fast version
# here we will:
#- Align
#- Record valid ligation events
#- sort pairsam file
#- remove PCR duplicates
#- generate pairs and bam files

# check if the file exists of it was created with a previous bam version 
outbam="${EDITED_DIR}/BAM/${filename}.PT.bam"
fileNotExistOrOlder "${outbam}" \
                    "${read1_path} ${read2_path}"

if [[ ${analyse} == "yes" ]]; then
    echo -e "Starting Alignment and duplicate removal ------------------------------\n"

    bwa mem -5SP -T0 -t${SLURM_CPUS_PER_TASK} ${REFERENCE_DIR}.fa ${read1_path} ${read2_path}| \
    \
    pairtools parse --min-mapq 40 --walks-policy 5unique \
    --max-inter-align-gap 30 --nproc-in ${SLURM_CPUS_PER_TASK} \
    --nproc-out ${SLURM_CPUS_PER_TASK} --chroms-path ${REFERENCE_DIR}.fa | \
    \
    pairtools sort --tmpdir=${tempdir} --nproc ${SLURM_CPUS_PER_TASK} | \
    \
    pairtools dedup --nproc-in ${SLURM_CPUS_PER_TASK} \
    --nproc-out ${SLURM_CPUS_PER_TASK} --mark-dups --output-stats "${statsOut}.txt" | \
    \
    pairtools split --nproc-in ${SLURM_CPUS_PER_TASK} \
    --nproc-out ${SLURM_CPUS_PER_TASK} --output-pairs ${outpair} --output-sam - | \
    \
    samtools view -bS -@${SLURM_CPUS_PER_TASK} | \
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

    # Get mapping stats
    # value thressholds for valid
    #Metric	                        Shallow Seq (20M)	Deep Seq (100-200M)
    #No-Dup Read Pairs	            >75%	            >50%
    #No-dup cis read pairs ≥ 1kb	>20%	            >20%

    python3 ${microCScripts}/00_NR_get_qc.py -p "${statsOut}.txt" > "${statsOut}_processed.txt"
    mv ${statsOut}_processed.txt ${statsOut}.txt

else
    echo -e "Alignment and duplicate removal - already done before--------------------------\n"

fi


###################################
# Chromosome filtering
#################################

# Remove chrM, chrUn, _random ... (remove sort.rmdup.rmblackls.bam)
echo -e "Starting remove chrM and useless chromosomes ------------------------------\n"


outbam2="${bamsPath}/${filename}.PT.rmchr.bam"

fileNotExistOrOlder "${outbam2}" \
                    "${outbam}"

if [[ ${analyse} == "yes" ]]; then
	samtools view -h ${outbam} | \
    awk '(!index($3, "random")) && (!index($3, "chrUn")) && ($3 != "chrM") && ($3 != "chrEBV")' | \
    samtools view -Sb - > ${outbam2}
    samtools index ${outbam2} -@ $SLURM_CPUS_PER_TASK
    echo -e "Remove chrM and useless chromosomes - done ----------------------------------\n"
    if [[ $removeTemp == 'yes' ]] ; then
        rm ${outbam}
    fi
else
    echo -e "Remove chrM and useless chromosomes - already done before--------------------------\n"
fi


###################################
# BigWigs 
#################################
# deeptools bamCoverage
# Normalize by CPM (This is the scaled bigWig we will use)

if [ ! -e ${basePath}/BigWig/ ]; then
	mkdir ${basePath}/BigWig/
fi

bigWigOut="${basePath}/BigWig/${filename}.PT.rmchr.bw"

fileNotExistOrOlder "${bigWigOut}" \
                    "${outbam2}"

if [[ ${analyse} == "yes" ]]; then
	bamCoverage --binSize 5 --normalizeUsing CPM --exactScaling \
    -b ${outbam2} -of bigwig \
    -o ${bigWigOut} --numberOfProcessors $SLURM_CPUS_PER_TASK
	
    echo -e "BigWig norm 1 - done ---------------------------------------------\n"

else
    echo -e "BigWig norm 1 - already done before ------------------------------\n"
fi



################################
# Preseq QC
################################

outPreseq="${basePath}/QC/preseq/${filename}.preseq"
fileNotExistOrOlder "${outPreseq}" \
                    "${outbam2}"

if [[ ${analyse} == "yes" ]]; then
	$preseq lc_extrap -bam -pe -extrap 2.1e9 -step 1e8 \
        -seg_len 1000000000 -output ${outPreseq} ${outbam2}

	
    echo -e "Preseq QC - done ---------------------------------------------\n"

else
    echo -e "Preseq QC - already done before ------------------------------\n"
fi



echo -e "END --------------------------------------------------"

seff $SLURM_JOBID

exit 0



# bamp="/home/jmendietaes/data/2021/microC/allProcessed/bamfiles"
# checkBams="DM-N_mC_6mn_1_S2.PT.rmchr.bam DM-F_mC_6mn_1_S4.PT.rmchr.bam \
# Mye_mC_3mn_1_S5.PT.rmchr.bam DM-N_mC_3mn_1_S1.PT.rmchr.bam Mye_mC_6mn_1_S6.PT.rmchr.bam"

# for bam in ${checkBams}; do
#     filename=$(echo $bam | sed 's/.PT.rmchr.bam//g')
#     outPreseq="${basePath}/QC/preseq/${filename}.preseq"
#     $preseq lc_extrap -bam -pe -extrap 2.1e9 -step 1e8 \
#             -seg_len 1000000000 -output ${outPreseq} ${bamp}/${bam}
# done