#!/bin/bash
# -*- ENCODING: UTF-8 -*-

##===============================================================================
## SLURM VARIABLES
#SBATCH --job-name=EMseq_fqToBw
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --time=02-00:00:00
#SBATCH -p medium
#SBATCH -o /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.out  
#SBATCH -e /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.err 
##SBATCH --dependency=afterany:793158


##SBATCH --mail-type=END
##SBATCH --mail-user=user@mail.es
# HOW TO RUN ME
# for i in *fastq.gz; do echo $i | sed 's/_R._001.fastq.gz//g' ; done | sort | uniq > samplesNames.txt  
# N=`cat samplesNames.txt | wc -l`
# sbatch --array=1-${N} /home/jmendietaes/programas/pipelines/EM-seq/cluster/01_EMseq_fqTObw.sh \
#/home/jmendietaes/data/2021/DNAme/sequencedData/NextSeq2000.RUN156.20230306 \
#/home/jmendietaes/data/2021/DNAme/allProcessed \
#/home/jmendietaes/referenceGenomes/mm10_reordered/mm10.reordered

# Nodes to avoid
# nodo11,nodo12,nodo06,nodo07 nodo08
# --nodelist
# The list may be specified as a comma-separated list of
#               hosts, a range of hosts (host[1-5,7,...]
# --exclude=<node name list>
# Explicitly exclude certain nodes from the  resources  granted  to  the
# job.

# module load BWA/0.7.17-foss-2018b
# /home/jmendietaes/programas/bwa-meth-master/bwameth.py --reference /home/jmendietaes/referenceGenomes/mm10_reordered/mm10.reordered.fa -p /home/jmendietaes/data/2021/DNAme/sequencedData/NextSeq2000.RUN156.20230306/pipelineOut/trimming/DM_DNAme1_1_S22.interleaved.fa -t 3 --read-group "@RG\\tID:CACTGTAG+AGTCGCTT\\tSM:DM_DNAme1_1_S22"

# SOURCE
# https://github.com/semenko/serpent-methylation-pipeline/tree/main
# https://github.com/nebiolabs/EM-seq/blob/master/em-seq.nf
# https://nf-co.re/methylseq/2.4.0/parameters#clip_r1

# UPDATES
# Array must be index 1 in my scenario
# it works with gziped fastq files
# not valid for atac seq, dosnt take into account the 9bp mapping sift (Tn5, hichikers paper)
# i can use nextflow or prepare time points in a file in case the code breaks
##===============================================================================
# Path to extra scripts
subScripts="/home/jmendietaes/programas/pipelines/EM-seq/cluster/sub-scripts"

## Trimming parameters
# Optimal trimming as defined by @nebiolabs and by checking Per base SC
# Also check https://github.com/FelixKrueger/Bismark/issues/509
trim_r1_5prime=8
trim_r1_3prime=5
# bs_seq has a huge adaptase bias
trim_r2_5prime=8
trim_r2_3prime=5
minimum_length=15
# Sequencing adapters, set to None to autoDetect 
adapter_r1="None"
adapter_r2="None"

## GLOBAL VARIABLES
PROJECT_DIR=$1
#PROJECT_DIR="/home/jmendietaes/data/2021/chip/NextSeq2000.RUN6.20210427"
# In my case, the fastq files are in project folder inside demux_fastq

RAW_FASTQ_DIR=$PROJECT_DIR"/demux_fastq"
# not used, but i will store here final files
#RESULTS_DIR=$PROJECT_DIR"/pipelineOut/results"
EDITED_DIR=$PROJECT_DIR"/pipelineOut"
FASTQ_DIR=$EDITED_DIR"/fastq" 
libMetricsP="${EDITED_DIR}/libMetrics"
# set to lowercase yes or not
removeTemp="yes"

## path to the main fully processed data output
basePath=$2
#basePath="/home/jmendietaes/data/2021/chip/allProcessed"
# here we will stored the final filtered bam files
bamsPath="${basePath}/bamfiles"

# Minimum number of nonconverted Cs on a read to consider it nonconverted
# https://github.com/nebiolabs/mark-nonconverted-reads
c_count=3

##  REFERENCE Genome info
REFERENCE_DIR=$3
#REFERENCE_DIR="/home/jmendietaes/referenceGenomes/mm10_reordered/mm10.reordered"
# REFERENCE_DIR includes the path and naming base of the reference genome, we will
#   add .sizes, .blacklist.bed to the base to get the rest of paths. Genome
#   indexes must follow the same base name

GenomeIndex=$REFERENCE_DIR
#chrOrder=$REFERENCE_DIR"/mm10_Bowtie2/names.txt"
chr_genome_size=$REFERENCE_DIR".sizes"
BlackList=$REFERENCE_DIR".blacklist.bed"
wigToBigWig="/home/jmendietaes/programas/binPath/wigToBigWig"
picardPath='/home/jmendietaes/programas/picard/picard.jar'
bwaMeth='/home/jmendietaes/programas/bwa-meth-master/bwameth.py'
fastpRun="/home/jmendietaes/programas/fastp/fastp"
MethylDackel="/home/jmendietaes/programas/MethylDackel/installation/MethylDackel"
#picardPath='$EBROOTPICARD/picard.jar'
##===============================================================================
## Required Software
#module load Trimmomatic/0.38-Java-1.8
#module load Trim_Galore/0.6.0
# trim_galore will requilre fastqc and cutadapt
# I have trim galore and cutadapt in conda
module load BWA/0.7.17-foss-2018b
module load SAMtools/1.12-GCC-10.2.0
module load Java/1.8.0_192
module load BEDTools/2.27.1-foss-2018b
module load HTSlib/1.11-GCC-10.2.0

#module load MACS2/2.1.0.20151222-foss-2018b-Python-2.7.15
#module load deepTools/3.2.0-foss-2018b-Python-2.7.15
# this is for bamCoverage, but i already have it in conda and newer version
#module load Homer/4.10-foss-2018b

##===============================================================================

echo julenmendieta92@gmail.com  jmendietaes@unav.es
echo Date ; date
echo script --basic-pre-process
echo LibraryType=PE
echo genome = ${GenomeIndex}
##===============================================================================

# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

##===============================================================================
if [ ! -e ${EDITED_DIR} ]; then
    mkdir -p ${EDITED_DIR}
fi

# no idea of what this is for
#if [[ ! -e ${EDITED_DIR} ]] && echo exists ; then
#    mkdir -p ${EDITED_DIR}
#fi

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
stepControl="${EDITED_DIR}/QC/pipelineStep_${filename}.txt"
summaryFile="${basePath}/QC/summary_${filename}.txt"
nCPU=$SLURM_CPUS_PER_TASK


if [[ ${adapter_r1} != "None" ]]; then 
    adapterStr="--adapter_sequence ${adapter_r1} --adapter_sequence_r2 ${adapter_r2}"
else
    adapterStr=""
fi

#############
# BEGIN: Create summary file
#############								
if [ ! -e ${EDITED_DIR}/QC/ ] && echo exists ; then
    mkdir -p ${EDITED_DIR}/QC/
fi

if [ ! -e ${basePath}/QC/ ] && echo exists ; then
    mkdir -p ${basePath}/QC/
fi

if [ ! -e ${stepControl} ] ; then
    touch ${stepControl}
fi

echo -e "Starting Summary file -------------------------------------- \n"

# check content of first line of step control file
linec=`sed "1q;d" ${stepControl}`
if [[ ${linec} != "Summary" ]]; then 
    echo -e "STARTING \n $(date) \n" 
    echo "SAMPLE: ${filename}" 
    

    # QC: read counts if file is gziped or not
    if file --mime-type ${read1_path} | grep -q gzip$; then
        Counts1="$(zcat ${read1_path} | echo $((`wc -l`/4)))"
        Counts2="$(zcat ${read2_path} | echo $((`wc -l`/4)))" 
    else
        echo "Not gzipped files"
        echo $read1_path
        #Counts1="$(cat ${RAW_FASTQ_DIR}/${filename}_R1_001.fastq | echo $((`wc -l`/4)))"
        #Counts2="$(cat ${RAW_FASTQ_DIR}/${filename}_R2_001.fastq| echo $((`wc -l`/4)))"
    fi


    echo -e "READ COUNTS" >> ${summaryFile}
    echo -e "sample name\tfastq name\tread count\tmillions" >> ${summaryFile}
    rc=$((Counts1/1000000))
    echo -e "${filename}\t${filename}_R1\t${Counts1}\t${rc}" >> ${summaryFile}
    rc=$((Counts2/1000000))
    echo -e "${filename}\t${filename}_R2\t${Counts2}\t${rc}\n" >> ${summaryFile}

    echo -e "Summary file - done -------------------------------------- \n"
    # store stage control info
    echo "Summary" > ${stepControl}
else
    echo -e "Summary file - already done before -------------------------------------- \n"
fi

#############
# TR: Trimming
#############
# We have to trim for repetitive seq, adapters, quality etc.
echo -e "Starting Trimming and FASTQC -------------------------------------- \n"

# check content of second line of step control file
linec=`sed "2q;d" ${stepControl}`
if [[ ${linec} != "Trim" ]]; then 
    # Trim_galore recommends less than 8 cpu
    # trimCPU=$([ $SLURM_CPUS_PER_TASK -le 7 ] && echo "$SLURM_CPUS_PER_TASK" \
    #         || echo "7")
    if [ ! -e ${EDITED_DIR}/fastQC/ ]; then
        mkdir -p ${EDITED_DIR}/trimming/  
        mkdir -p ${EDITED_DIR}/fastQC/
    fi

    if [ ! -e ${EDITED_DIR}/fastQC/${filename}_R2_001.fastq_fastqc.html ]; then
        #trim_galore --cores $trimCPU --paired --fastqc --gzip \
        #--output_dir ${EDITED_DIR}/trimming/ --fastqc_args "--outdir ${EDITED_DIR}/fastQC/" \
        #${read1_path} ${read2_path}


        # ${fastpRun} --stdin --stdout -l 2 -Q ${trim_polyg} --interleaved_in --overrepresentation_analysis \
        #     -j "!{fq_set.library}_fastp.json" 2> fastp.stderr

        # inst_name=$(zcat -f ${read1_path} | head -n 1 | cut -f 1 -d ':' | sed 's/^@//')
        # fastq_barcode=$(zcat -f ${read1_path} | head -n 1 | sed -r 's/.*://')

        # if [[ "${inst_name:0:2}" == 'A0' ]] || [[ "${inst_name:0:2}" == 'NS' ]] || \
        #     [[ "${inst_name:0:2}" == 'NB' ]] || [[ "${inst_name:0:2}" == 'VH' ]] ; then
        #     trim_polyg='--trim_poly_g'
        #     echo '2-color instrument: poly-g trim mode on'
        # else
        #     trim_polyg=''
        # fi


        # trim_poly_g is enabled for NextSeq/NovaSeq data by default
        ${fastpRun} --in1 ${read1_path} --in2 ${read2_path} \
        --trim_front1 ${trim_r1_5prime} --trim_tail1 ${trim_r1_3prime} \
        --trim_front2 ${trim_r2_5prime} --trim_tail2 ${trim_r2_3prime} \
        --length_required ${minimum_length} ${adapterStr} \
        --json ${EDITED_DIR}/fastQC/${filename}.fastp.json \
        --html ${EDITED_DIR}/fastQC/${filename}.fastp.html \
        --thread ${nCPU} --verbose \
        --failed_out ${EDITED_DIR}/fastQC/${filename}.fastp.failed.fa.gz \
        --stdout 2>${EDITED_DIR}/fastQC/${filename}.fastp.out \
        > ${EDITED_DIR}/trimming/${filename}.interleaved.fa

        readsPass=$(grep "reads passed filter" ${EDITED_DIR}/fastQC/${filename}.fastp.out)
        readsPass=$(echo ${readsPass} | sed 's/reads passed filter: //g')
        readsPass=$(($readsPass/2))
        echo -e "FASTP TRIMMING" >> ${summaryFile}
        echo -e "reads passing filters: ${readsPass}" >> ${summaryFile}

    # grep "reads passed filter" ${EDITED_DIR}/fastQC/${filename}.fastp.out
    # dividir entre 2

    fi
    echo -e "Trimming - done -------------------------------------- \n"

    # store stage control info
    echo "Trim" >> ${stepControl}

else
    echo -e "Trimming - already done before -------------------------------------- \n"
fi

# define new path variables
trimedReads="${EDITED_DIR}/trimming/${filename}.interleaved.fa"

#############
#Aling # Mapping
############
echo -e "Starting Alignment -------------------------------------- \n"

if [ ! -e ${EDITED_DIR}/BAM/ ]; then
    mkdir -p ${EDITED_DIR}/BAM/
    mkdir -p ${EDITED_DIR}/unMapped/
fi

samFile="${EDITED_DIR}/BAM/${filename}.sam"

# check content of third line of step control file
extr_unmap="${EDITED_DIR}/unMapped/${filename}_R%_001_unmap.fastq.gz"
linec=`sed "3q;d" ${stepControl}`
if [[ ${linec} != "Align" ]]; then 

    # ID: In Illumina data, read group IDs are composed using the flowcell name 
    #    and lane number, making them a globally unique identifier across all 
    #    sequencing data in the world
    # SM: The name of the sample sequenced in this read group
    fastq_barcode=$(zcat -f ${read1_path} | head -n 1 | sed -r 's/.*://')
    #zcat -f  ${read1_path} | head -n 1 | cut -f 1 -d ':' | sed 's/^@//'
    # Last time that worked i tried with "@RG\\tID:${filename}\\tSM:${filename}"
    read_group="@RG\\tID:${fastq_barcode}\\tSM:${filename}"
    
    # CHECK: 
    # This only works with one CPU in my cluster, "Broken Pipe" error
    # Some people change line 362 of bwameth.py to inlcude "-1" option in bwa mem
    # this hidden option -1 still enables multi-threading for mapping, but 
    # disables fastq reading and sam writing on different threads. 
    # Didn't work for me :(
    # https://github.com/lh3/bwa/issues/102

    ${bwaMeth} -p -t 1 --read-group ${read_group} \
        --reference ${GenomeIndex}.fa ${trimedReads} \
        2> ${EDITED_DIR}/BAM/${filename}.bwameth_err \
        > ${samFile}

        
    echo -e "Alignment - done -------------------------------------- \n"
    # store stage control info
    echo "Align" >> ${stepControl}

    # remove temporal trimmed fastq files
    if [[ $removeTemp == 'yes' ]] ; then
        rm ${trimedReads}
    fi
else
    echo -e "Alignment - already done before ------------------------------ \n"
fi

# SAM to BAM (remove SAM)
echo -e "Starting SAM to BAM -------------------------------------- \n"

bamFile="${EDITED_DIR}/BAM/${filename}.bam"
bamSort="${EDITED_DIR}/BAM/${filename}.sort.bam"
# check content of fourth line of step control file
linec=`sed "4q;d" ${stepControl}`
if [[ ${linec} != "SamBam" ]]; then 
    samtools view -o ${bamFile} -bhS -@ ${nCPU} \
                    ${samFile} 
    samtools sort -o ${bamSort} ${bamFile} 
    samtools index ${bamSort} -@ ${nCPU}

    echo -e "SAM to BAM - done -------------------------------------- \n"
    # store stage control info
    echo "SamBam" >> ${stepControl}

    # delete intermediate files
    if [[ $removeTemp == 'yes' ]] ; then
        rm ${samFile}
        rm ${bamFile}
    fi

else
    echo -e "SAM to BAM - already done before ------------------------------ \n"
fi




##############################
#Dupicates marking and removal
##############################
echo -e "Starting mark and remove duplicates ------------------------------------------------ \n"
# Samtools flags to remove duplicates + unmapped reads:
#   a) -F 1028 => read unmapped // duplicate
#   b) -F 1804 => read unmapped // mate unmapped // not primary alignment // read failing platform // duplicate

if [ ! -e ${libMetricsP} ]; then
    mkdir -p ${libMetricsP}
fi

bamSortMarkDup="${EDITED_DIR}/BAM/${filename}.sort.markdup.bam"
bamSortRmDup="${EDITED_DIR}/BAM/${filename}.sort.rmdup.bam"

# check content of sixth line of step control file
linec=`sed "5q;d" ${stepControl}`
if [[ ${linec} != "markDup" ]]; then 
   # Identify duplicated reads by flaging them in a new ioutput BAM file
    java -Xmx24G -jar ${picardPath} MarkDuplicates I=${bamSort} O=${bamSortMarkDup} \
        R=${GenomeIndex}.fa \
        METRICS_FILE=${libMetricsP}/${filename}.bam.mkdup.txt
    # remove reads marked as duplicated
    samtools view -o ${bamSortRmDup} -@ ${nCPU} -bh -F 1804 ${bamSortMarkDup}
    
    # remove orphan reads?
    # bampe_rm_orphan.py ${bam[0]} ${prefix}.bam --only_fr_pairs

    # QC: Show duplicates in samtools flags
    echo -e "\nSAMTOOLS FLAGSTAT - DUPLICATES" >> ${summaryFile}

    samtools flagstat ${bamSortMarkDup} >> ${summaryFile}

    echo -e "\n"  >> ${summaryFile}
   
    if [[ $removeTemp == 'yes' ]] ; then
        rm ${bamSort}
    fi
    echo -e "Remove duplicates - done ------------------------------------------------ \n"
    # store stage control info
    echo "markDup" >> ${stepControl}
else
    echo -e "Remove duplicates - already done before ------------------------------------ \n"

fi 

echo -e "DATE CHECK \n $(date) \n" 



##############################
# Remove chrM, chrUn, _random ... (remove sort.rmdup.rmblackls.bam)
##############################

echo -e "Starting remove chrM and useless chromosomes ------------------------------\n"

if [ ! -e ${bamsPath} ]; then
    mkdir -p ${bamsPath}/allRead
    mkdir -p ${bamsPath}/onlyConverted
fi

bamSortRmDupBlackChr="${bamsPath}/allRead/${filename}.sort.rmdup.rmchr.bam"

# check content of eigth line of step control file
linec=`sed "6q;d" ${stepControl}`
if [[ ${linec} != "Remove" ]]; then 
	samtools view -h ${bamSortRmDup} | \
    awk '(!index($3, "random")) && (!index($3, "chrUn")) && ($3 != "chrM") && ($3 != "chrEBV")' | \
    samtools view -Sb - > ${bamSortRmDupBlackChr}
    samtools index ${bamSortRmDupBlackChr} -@ ${nCPU}

    # QC: Show final reads
    echo -e "\nSAMTOOLS FLAGSTAT - FINAL READS" >> ${summaryFile}
    samtools flagstat ${bamSortRmDupBlackChr} >> ${summaryFile}
    echo -e "\n" >> ${summaryFile}

    echo -e "Remove chrM and useless chromosomes - done ----------------------------------\n"
    # store stage control info
    echo "Remove" >> ${stepControl}
    # if [[ $removeTemp == 'yes' ]] ; then
    #     rm ${bamSortMarkDupBlack}
    # fi
else
    echo -e "Remove chrM and useless chromosomes - already done before--------------------------\n"
fi

# awk '(!index($1, "random")) && (!index($1, "chrUn")) && ($1 != "chrM") && ($1 != "chrEBV")' mm10.reordered.sizes


###################################
# Mark unconverted reads
#################################

echo -e "Starting mark of unconverted reads ------------------------------\n"

if [ ! -e ${bamsPath} ]; then
    mkdir -p ${bamsPath}
fi

samSortRmDupBlackChrConv="${EDITED_DIR}/BAM/${filename}.sort.rmdup.rmchr.convrtd.sam"
bamSortRmDupBlackChrConv="${EDITED_DIR}/BAM/${filename}.sort.rmdup.rmchr.convrtd.bam"
# check content of eigth line of step control file
linec=`sed "7q;d" ${stepControl}`
if [[ ${linec} != "convrtd" ]]; then 
    
    ${subScripts}/mark-nonconverted-reads.py --reference ${GenomeIndex}.fa \
                --bam ${bamSortRmDupBlackChr} --out ${samSortRmDupBlackChrConv} \
                --c_count ${c_count} --flag_reads 2> ${EDITED_DIR}/BAM/log.nonconverted.txt
    
    samtools view -o ${bamSortRmDupBlackChrConv} -bhS -@ ${nCPU} \
                    ${samSortRmDupBlackChrConv} 

    echo -e "mark of unconverted reads - done ----------------------------------\n"
    # store stage control info
    echo "convrtd" >> ${stepControl}

else
    echo -e "mark of unconverted reads - already done before--------------------------\n"
fi

#samtools view -F 512 ${samSortRmDupBlackChrConv}
#samtools view -o ${bamSortRmDup} -@ $SLURM_CPUS_PER_TASK -bh -F 512 ${bamSortMarkDup}


###################################
# Filter out unconverted reads
#################################
bamRmDupBlackChrConvFlt="${EDITED_DIR}/BAM/${filename}.rmdup.rmchr.convrtdFlt.bam"
bamSortRmDupBlackChrConvFlt="${bamsPath}/onlyConverted/${filename}.sort.rmdup.rmchr.convrtdFlt.bam"
echo -e "Starting filtering of unconverted reads ------------------------------\n"
# check content of eigth line of step control file
linec=`sed "8q;d" ${stepControl}`
if [[ ${linec} != "convrtdflt" ]]; then 
    
    samtools view -o ${bamRmDupBlackChrConvFlt} -@ $SLURM_CPUS_PER_TASK \
                    -bh -F 512 ${bamSortRmDupBlackChrConv}
    samtools sort -o ${bamSortRmDupBlackChrConvFlt} ${bamRmDupBlackChrConvFlt} 
    samtools index ${bamSortRmDupBlackChrConvFlt} -@ ${nCPU}

    # Count number of reads
    echo -e "\nSAMTOOLS FLAGSTAT - FINAL CONVERTED READS" >> ${summaryFile}
    samtools flagstat ${bamSortRmDupBlackChrConvFlt} >> ${summaryFile}
    echo -e "\n" >> ${summaryFile}

    echo -e "filtering of unconverted reads - done ----------------------------------\n"
    # store stage control info
    echo "convrtdflt" >> ${stepControl}

else
    echo -e "filtering of unconverted reads - already done before--------------------------\n"
fi


###################################
# BigWigs 
#################################
# deeptools bamCoverage
# Normalize by CPM (This is the scaled bigWig we will use)
echo -e "Starting BigWigs --------------------------------------------------\n"

if [ ! -e ${basePath}/BigWig/allRead ]; then
	mkdir ${basePath}/BigWig/allRead
    mkdir ${basePath}/BigWig/onlyConverted
fi

# if [ ! -e ${EDITED_DIR}/BigWig/ ]; then
# 	mkdir ${EDITED_DIR}/BigWig/
# fi

bigWigOutAll="${basePath}/BigWig/allRead/${filename}.sort.rmdup.rmchr.norm.bw"
bigWigOutCvrtd="${basePath}/BigWig/onlyConverted/${filename}.sort.rmdup.rmchr.convrtdFlt.norm.bw"

# check content of tenth line of step control file
linec=`sed "9q;d" ${stepControl}`
if [[ ${linec} != "BigWnorm1" ]]; then 
	bamCoverage --binSize 5 --normalizeUsing CPM --exactScaling \
    -b ${bamSortRmDupBlackChr} -of bigwig \
    -o ${bigWigOutAll} --numberOfProcessors ${nCPU}


    bamCoverage --binSize 5 --normalizeUsing CPM --exactScaling \
    -b ${bamSortRmDupBlackChrConvFlt} -of bigwig \
    -o ${bigWigOutCvrtd} --numberOfProcessors ${nCPU}
	# bamCoverage --binSize 20 --normalizeUsing RPKM --effectiveGenomeSize 2913022398 \
	# -b ${BaseFolder}BAM/${filename}.sort.rmdup.rmblackls.rmchr.bam -of bigwig \
	# -o ${BaseFolder}BigWig/${filename}.sort.rmdup.rmblackls.rmchr.norm.bw
    
    echo -e "BigWig norm 1 - done ---------------------------------------------\n"
    # store stage control info
    echo "BigWnorm1" >> ${stepControl}
else
    echo -e "BigWig norm 1 - already done before ------------------------------\n"
fi



# I might remove this from the pipeline in the future
# echo -e "Starting BigWigs different normalization---------------------------\n"

# bigWigOut2="${EDITED_DIR}/BigWig/${filename}.sort.rmdup.rmblackls.rmchr.gencovnorm.bw"

# # check content of eleventh line of step control file
# linec=`sed "11q;d" ${stepControl}`
# if [[ ${linec} != "BigWnorm2" ]]; then 
#     samtools view -b ${bamSortMarkDupBlackChr} | \
#             genomeCoverageBed -ibam stdin -g $chr_genome_size -bg | \
#             $wigToBigWig -clip stdin $chr_genome_size ${bigWigOut2}
#     #samtools view -b ${bamSortMarkDupBlackChr} | genomeCoverageBed -ibam stdin -g $chr_genome_size -bg > ${basePath}/BigWig/${filename}.sort.rmdup.rmblackls.rmchr.gencovnorm.bedGraph
 
#     # QC: print scaling factors
#     echo -e "BIGWIG SCALING FACTOR for first: CPM binSize 5" >> ${summaryFile}
#     echo -e "BIGWIG SCALING FACTOR for second normalization: Calculated with genomeCoverageBed " >> ${summaryFile}
#     echo -e "\n" >> ${summaryFile}

#     echo -e "BigWigs second - done ---------------------------------------------\n"
#     # store stage control info
#     echo "BigWnorm2" >> ${stepControl}
# else
#     echo -e "BigWigs second - already done before ------------------------------\n"
# fi    
#sortBed -i "$i" | genomeCoverageBed -bg -i stdin -g ~/alignments/chr_genome_size/chr_hg38_size/hg38.chrom.sizes | wigToBigWig -clip stdin ~/alignments/chr_genome_size/chr_hg38_size/hg38.chrom.sizes ${f%.*PE*}.bw


###################################
# MethylDackel extract 
#################################

if [ ! -e ${basePath}/methylDackel/ ]; then
	mkdir ${basePath}/methylDackel/
fi

cd ${basePath}/methylDackel/
echo -e "Starting MethylDackel extract ---------------------------------------\n"
# check content of tenth line of step control file
linec=`sed "10q;d" ${stepControl}`
if [[ ${linec} != "methylExtract" ]]; then 


    ${MethylDackel} extract --methylKit -@ ${nCPU} \
            --CHH --CHG -o ${basePath}/methylDackel/${filename} ${GenomeIndex}.fa ${bamSortRmDupBlackChr}
    ${MethylDackel} extract -@ ${nCPU} \
            --CHH --CHG -o ${basePath}/methylDackel/${filename} ${GenomeIndex}.fa ${bamSortRmDupBlackChr}
    pigz -p ${nCPU} ${filename}*.methylKit
    pigz -p ${nCPU} ${filename}*.bedGraph

    echo -e "MethylDackel extract - done ---------------------------------------------\n"
    # store stage control info
    echo "methylExtract" >> ${stepControl}
else
    echo -e "MethylDackel extract - already done before ------------------------------\n"
fi


###################################
# MethylDackel mbias 
#################################

echo -e "Starting MethylDackel mbias ---------------------------------------\n"
# check content of tenth line of step control file
linec=`sed "10q;d" ${stepControl}`
if [[ ${linec} != "methylmbias" ]]; then 

    ${MethylDackel} mbias -@ ${nCPU} ${GenomeIndex}.fa ${bamSortRmDupBlackChr} \
        ${basePath}/methylDackel/${filename} --txt \
        > ${basePath}/methylDackel/${filename}.mbias.txt


    echo -e "MethylDackel mbias - done ---------------------------------------------\n"
    # store stage control info
    echo "methylmbias" >> ${stepControl}
else
    echo -e "MethylDackel mbias - already done before ------------------------------\n"
fi



# Continue analysis with
#https://bioconductor.org/packages/release/bioc/vignettes/methylKit/inst/doc/methylKit.html

echo -e "END --------------------------------------------------"

seff $SLURM_JOBID

exit 0


# A few approaches can mitigate the observability issue, including filtering 
# for positions that have at least a minimal depth (e.g., 3×, 5×, etc.) 
# and using algorithms that account for read depth when comparing cytosines 
# between samples.

# Because methylation data are inherently bimodal (i.e., most β scores are near 
# 0 or 1, as explored in Figures 4A–4D), methods that use the binomial or 
# β-binomial distribution tend to exhibit better performance for methylation 
# data than statistical tests that use other distributions