#!/bin/bash
# -*- ENCODING: UTF-8 -*-

##===============================================================================
## SLURM VARIABLES
#SBATCH --job-name=cutTag_fqToBw
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=00-08:00:00
#SBATCH -p short
#SBATCH -o /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.out  
#SBATCH -e /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.err 
##SBATCH --dependency=afterany:805719

##SBATCH --mail-type=END
##SBATCH --mail-user=user@mail.es
# HOW TO RUN ME
# for i in *fastq.gz; do echo $i | sed 's/_R._001.fastq.gz//g' ; done | sort | uniq > samplesNames.txt  
# N=`cat samplesNames.txt | wc -l`
# sbatch --array=1-${N} /home/jmendietaes/programas/pipelines/cutTag/cluster/01_cutTag_fqToBw.sh \
#/home/jmendietaes/data/2021/cutTag/sequencedData/NextSeq2000.RUN108.221005 \
#/home/jmendietaes/data/2021/cutTag/allProcessed \
#/home/jmendietaes/referenceGenomes/mm10_reordered/mm10.reordered

#DM-Kmt2a_H3K4me1_CUTTAG_2_S9

# PURPOSE
# Analyse Cut and Tag data generated with Spike-in sequences like the ones from 
# SNAP-CUTANA K-MetStat

# Path to base name to point to .fa and .EpiListRenamed
# .fa would be the fasta file with the spike-in sequences
# .EpiListRenamed contains a list of in-house ids for epigenetic marks
#   that are included at CUTANA_K-MetStat
spikeInRef="/home/jmendietaes/referenceGenomes/spikeIn/SNAP-CUTANA_K-MetStat"

# Control ID (to remove duplicates if its IgG for example)
controlID=IgG
# UPDATES
# Code adapted from https://yezhengstat.github.io/CUTTag_tutorial/ and
# https://github.com/nf-core/cutandrun taking into account we donw use drosophila
# spike-In but SNAP-CUTANA K-MetStat
# https://www.epicypher.com/products/nucleosomes/snap-cutana-k-metstat-panel

##===============================================================================
## GLOBAL VARIABLES
PROJECT_DIR=$1
#PROJECT_DIR="/home/jmendietaes/data/2021/cutTag/sequencedData/NextSeq2000.RUN108.221005"
# In my case, the fastq files are in project folder inside demux_fastq

RAW_FASTQ_DIR=$PROJECT_DIR"/demux_fastq"
# not used, but i will store here final files
#RESULTS_DIR=$PROJECT_DIR"/pipelineOut/results"
EDITED_DIR=$PROJECT_DIR"/pipelineOut"
FASTQ_DIR=$EDITED_DIR"/fastq" 
libMetricsP="${EDITED_DIR}/libMetrics"
# set to lowercase yes or not
removeTemp="no"

## path to the main fully processed data output
basePath=$2
#basePath="/home/jmendietaes/data/2021/cutTag/allProcessed"
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
chr_genome_size=$REFERENCE_DIR".sizes"
BlackList=$REFERENCE_DIR".blacklist.bed"
wigToBigWig="/home/jmendietaes/programas/binPath/wigToBigWig"
picardPath='/home/jmendietaes/programas/picard/picard.jar'
# Path to Qualimap
qualimap="/home/jmendietaes/programas/pipelines/qualimap/qualimap_v2.2.1/qualimap"
#picardPath='$EBROOTPICARD/picard.jar'
##===============================================================================
## Required Software
#module load Trimmomatic/0.38-Java-1.8
#module load Trim_Galore/0.6.0
# trim_galore will requilre fastqc and cutadapt
# I have trim galore and cutadapt in conda
module load FastQC/0.11.8-Java-1.8
module load Bowtie2/2.3.4.2-foss-2018b
module load SAMtools/1.12-GCC-10.2.0
module load picard/2.18.17-Java-1.8
module load R/4.0.0-foss-2018b
module load Java/1.8.0_192
module load BEDTools/2.27.1-foss-2018b
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


# Load assigned RAM
adjustedMem=$(echo "print(int(round($SLURM_MEM_PER_NODE*0.96,0)/1000))" | python3)
# fast check in cases where RAM given in Gb (value will be lower than 1 in most cases)
# 100Gb / 1000 = 0.1
if [ "$adjustedMem" -lt "5" ]; then
    adjustedMem=$(echo "print(int(round($SLURM_MEM_PER_NODE*0.96,0)))" | python3)
fi  
# Convert to bytes
adjustedMem_bytes=$(printf ${adjustedMem} | numfmt --from-unit=Mi)
adjustedMem_Gb=${adjustedMem}
echo  "Memory in Gb:"
echo ${adjustedMem}

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

epiName=(${filename//_/ }); epiName=${epiName[1]}; 
epiName=(${epiName//-/ }); epiName=${epiName[0]}


#############
# BEGIN: Create summary file
#############								
if [ ! -e ${EDITED_DIR}/QC/ ] && echo exists ; then
    mkdir -p ${EDITED_DIR}/QC
    mkdir -p ${EDITED_DIR}/FASTQ/
fi

if [ ! -e ${basePath}/QC/qualimap/${filename}/ ] && echo exists ; then
    mkdir -p ${basePath}/QC/qualimap/${filename}/
fi

if [ ! -e ${stepControl} ] ; then
    touch ${stepControl}
fi

bamSTR=""
echo -e "Starting Summary file -------------------------------------- \n"

# check content of first line of step control file
linec=`sed "1q;d" ${stepControl}`
if [[ ${linec} != "Summary" ]]; then 
    echo -e "STARTING \n $(date) \n" 

    # Is more efficient to keep and unzipped copy of the fastq files for the 
    # Spike in check
    read1_nogz="${EDITED_DIR}/FASTQ/${filename}_R1_001.fastq"
    read2_nogz="${EDITED_DIR}/FASTQ/${filename}_R2_001.fastq"

    pigz --decompress --keep -p $SLURM_CPUS_PER_TASK \
            --stdout ${read1_path} > ${read1_nogz}
    pigz --decompress --keep -p $SLURM_CPUS_PER_TASK \
            --stdout ${read2_path} > ${read2_nogz}

    echo "SAMPLE: ${filename}" 
    

    # # QC: read counts if file is gziped or not
    # if file --mime-type ${read1_path} | grep -q gzip$; then
    #     Counts1="$(zcat ${read1_path} | echo $((`wc -l`/4)))"
    #     Counts2="$(zcat ${read2_path} | echo $((`wc -l`/4)))" 
    # else
    #     echo "Not gzipped files"
    #     echo $read1_path
    Counts1="$(cat ${read1_nogz} | echo $((`wc -l`/4)))"
    Counts2="$(cat ${read2_nogz} | echo $((`wc -l`/4)))"
    # fi


    echo -e "READ COUNTS" >> ${summaryFile}
    echo -e "sample name\tfastq name\tread count\tmillions" >> ${summaryFile}
    rc=$((Counts1/1000000))
    echo -e "${filename} \t ${filename}_R1 \t ${Counts1} \t ${rc}" >> ${summaryFile}
    rc=$((Counts2/1000000))
    echo -e "${filename} \t ${filename}_R2 \t ${Counts2} \t ${rc} \n" >> ${summaryFile}

    echo -e "Summary file - done -------------------------------------- \n"
    # store stage control info
    echo "Summary" > ${stepControl}
else
    echo -e "Summary file - already done before -------------------------------------- \n"
fi

#############
# TR: Trimming - Recommended to be skipped in <25bp read sequencing
#############
# We have to trim for repetitive seq, adapters, quality etc.
echo -e "Starting Trimming and FASTQC -------------------------------------- \n"

# check content of second line of step control file
linec=`sed "2q;d" ${stepControl}`
if [[ ${linec} != "Trim" ]]; then 
    # Trim_galore recommends less than 8 cpu
    trimCPU=$([ $SLURM_CPUS_PER_TASK -le 7 ] && echo "$SLURM_CPUS_PER_TASK" \
                || echo "7")
    if [ ! -e ${EDITED_DIR}/fastQC/ ]; then
        mkdir -p ${EDITED_DIR}/trimming/  
        mkdir -p ${EDITED_DIR}/fastQC/
    fi

    if [ ! -e ${EDITED_DIR}/fastQC/${filename}_R2_001.fastq_fastqc.html ]; then
        trim_galore --cores $trimCPU --paired --fastqc --gzip \
        --output_dir ${EDITED_DIR}/trimming/ --fastqc_args "--outdir ${EDITED_DIR}/fastQC/" \
        ${read1_path} ${read2_path}
    fi
    echo -e "Trimming - done -------------------------------------- \n"

    # store stage control info
    echo "Trim" >> ${stepControl}

else
    echo -e "Trimming - already done before -------------------------------------- \n"
fi

# define new path variables
trimedRead1="${EDITED_DIR}/trimming/${filename}_R1_001_val_1.fq.gz"
trimedRead2="${EDITED_DIR}/trimming/${filename}_R2_001_val_2.fq.gz"

#############
# Aling to reference genome and Spike-in
############
echo -e "Starting Alignment -------------------------------------- \n"

if [ ! -e ${EDITED_DIR}/BAM/ ]; then
    mkdir -p ${EDITED_DIR}/BAM/
    mkdir -p ${EDITED_DIR}/unMapped/
fi

if [ ! -e ${basePath}/spikeQC/ ] && echo exists ; then
    mkdir -p ${basePath}/spikeQC/counts
fi

samFile="${EDITED_DIR}/BAM/${filename}.sam"

# check content of third line of step control file
extr_unmap="${EDITED_DIR}/unMapped/${filename}_R%_001_unmap.fastq.gz"
linec=`sed "3q;d" ${stepControl}`
if [[ ${linec} != "Align" ]]; then 
    ## The spike-in sequences
    barcodes=$(grep -v ">" $spikeInRef.fa)
    # bowtie2 -p $SLURM_CPUS_PER_TASK \
    #     --end-to-end --very-sensitive --phred33 \
    #     -x ${spikeInRef} -1 ${read1_path} -2 ${read2_path} \
    #     -S ${EDITED_DIR}/BAM/${filename}_spikeIn.sam
    # According to EpiCypher their sequences are searched with a grep

    if [ $SLURM_CPUS_PER_TASK -ge 3 ] ; then
        for bc in ${barcodes} ; do
            grep -c $bc ${read1_nogz}
        done > ${basePath}/spikeQC/counts/${filename}_R1.txt &

        for bc in ${barcodes} ; do
            grep -c $bc ${read2_nogz}
        done > ${basePath}/spikeQC/counts/${filename}_R2.txt &

        # Remove two cpus for spike check, but keep at least 1
        useCPU=$(($SLURM_CPUS_PER_TASK - 2))
        useCPU=$([ $useCPU -ge 1 ] && echo "$useCPU" \
                || echo "1")
    else
        for bc in ${barcodes} ; do
            grep -c $bc ${read1_nogz}
        done > ${basePath}/spikeQC/counts/${filename}_R1.txt

        for bc in ${barcodes} ; do
            grep -c $bc ${read2_nogz}
        done > ${basePath}/spikeQC/counts/${filename}_R2.txt

        useCPU=$SLURM_CPUS_PER_TASK
    fi
    
    ## To the main genome 
    bowtie2 -p $useCPU \
    --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 \
    -x $GenomeIndex -1 ${trimedRead1} -2 ${trimedRead2} \
    -S ${samFile} >> ${summaryFile} 2>&1
        
    
    # Make sure alignment of spike in ended before exiting this stage
    wait $(jobs -rp)

    echo -e "Alignment - done -------------------------------------- \n"
    # store stage control info
    echo "Align" >> ${stepControl}

    # remove temporal trimmed fastq files
    if [[ $removeTemp == 'yes' ]] ; then
        rm ${trimedRead1}
        rm ${trimedRead2}
    fi
else
    echo -e "Alignment - already done before ------------------------------ \n"
fi


# SAM to BAM (remove SAM)
echo -e "Starting SAM to BAM -------------------------------------- \n"

bamFile="${EDITED_DIR}/BAM/${filename}.bam"

# check content of fourth line of step control file
linec=`sed "4q;d" ${stepControl}`
if [[ ${linec} != "SamBam" ]]; then 
    samtools view -o ${bamFile} -bhS -@ $SLURM_CPUS_PER_TASK \
                    ${samFile} 

    echo -e "SAM to BAM - done -------------------------------------- \n"
    # store stage control info
    echo "SamBam" >> ${stepControl}
else
    echo -e "SAM to BAM - already done before ------------------------------ \n"
fi


# Sort BAM by location (remove original BAM)

echo -e "Starting BAM sorting --------------------------------------\n"

bamSTR="sort"
bamSort="${EDITED_DIR}/BAM/${filename}.${bamSTR}.bam"

# check content of fifth line of step control file
linec=`sed "5q;d" ${stepControl}`
if [[ ${linec} != "BamSort" ]]; then 
    samtools sort -@ $SLURM_CPUS_PER_TASK -o ${bamSort} ${bamFile} 
    echo -e "BAM sorting - done --------------------------------------\n"
    # store stage control info
    echo "BamSort" >> ${stepControl}

    # delete intermediate files
    if [[ $removeTemp == 'yes' ]] ; then
        rm ${samFile}
        rm ${bamFile}
    fi
else
    echo -e "BAM sorting - already done before -----------------------------\n"
fi


##############################
# Dupicates marking (Removal not recommended in cut&Tag unless IgG)
##############################
# Control samples such as those from IgG datasets have relatively high 
#   duplication rates due to non-specific interactions with the genome; 
#   therefore, it is appropriate to remove duplicates from control samples.

echo -e "Starting mark duplicates ------------------------------------------------ \n"
# Samtools flags to remove duplicates + unmapped reads:
#   a) -F 1028 => read unmapped // duplicate
#   b) -F 1804 => read unmapped // mate unmapped // not primary alignment 
#// read failing platform // duplicate
#   c) -F 780 =>  read unmapped // mate unmapped // not primary alignment 
#// read failing platform 

if [ ! -e ${libMetricsP} ]; then
    mkdir -p ${libMetricsP}
fi


bamSortMarkDup="${EDITED_DIR}/BAM/${filename}.sort.markdup.bam"
bamSortRmUnM="${EDITED_DIR}/BAM/${filename}.sort.rmUnM.bam"

# Check if this is a control that needs to be de-duplicated
filter="no"
if [[ $filename == *_${controlID}_* ]] || [[ $filename == *_${controlID}-* ]]; then
    filter="yes"
    bamSTR="${bamSTR}.rmdup"
    bamSortFlt="${EDITED_DIR}/BAM/${filename}.${bamSTR}.bam"
else
    bamSTR="${bamSTR}.rmUnM"
    bamSortFlt=${bamSortRmUnM}
fi



# check content of sixth line of step control file
linec=`sed "6q;d" ${stepControl}`
if [[ ${linec} != "DupCheck" ]]; then 
    # Estimate the numbers of unique molecules in a sequencing library
    java -Xmx24G -jar ${picardPath} EstimateLibraryComplexity I=${bamSort} \
        O=${libMetricsP}/${filename}.bam.lib_comp.txt \
        VALIDATION_STRINGENCY=SILENT
    # Check insert size distribution and read orientation
    java -Xmx24G -jar ${picardPath} CollectInsertSizeMetrics I=${bamSort} \
        O=${libMetricsP}/${filename}.bam.insert_size.txt \
        H=${libMetricsP}/${filename}.bam.insert_size.pdf \
        VALIDATION_STRINGENCY=SILENT
    # Identify duplicated reads by flaging them in a new ioutput BAM file
    java -Xmx24G -jar ${picardPath} MarkDuplicates I=${bamSort} \
        O=${bamSortMarkDup} \
        METRICS_FILE=${libMetricsP}/${filename}.bam.mkdup.txt

    if [[ $filter == "yes" ]]; then
        # remove reads marked as duplicated and non-aligned
        samtools view -o ${bamSortFlt} -@ $SLURM_CPUS_PER_TASK -bh \
                    -F 1804 ${bamSortMarkDup}
    else
        # remove reads marked as non-aligned
        samtools view -o ${bamSortFlt} -@ $SLURM_CPUS_PER_TASK -bh \
                    -F 780 ${bamSortMarkDup}
    fi
    
    # Check insert size distribution and read orientation after removing duplicates
    java -Xmx24G -jar ${picardPath} CollectInsertSizeMetrics I=${bamSortFlt} \
        O=${libMetricsP}/${filename}.rmdup.bam.insert_size.txt \
        H=${libMetricsP}/${filename}.rmdup.bam.insert_size.pdf \
        VALIDATION_STRINGENCY=SILENT

    # remove orphan reads?
    # bampe_rm_orphan.py ${bam[0]} ${prefix}.bam --only_fr_pairs

    # QC: Show duplicates in samtools flags
    echo -e "\nSAMTOOLS FLAGSTAT - DUPLICATES" >> ${summaryFile}
    echo -e "~ >70% uniquely mapped reads expected\nWith ~ <50% \
be cautious" >> ${summaryFile}

    samtools flagstat ${bamSortMarkDup} >> ${summaryFile}

    # Add % of duplicates from aligned reads
    mapped=$(grep -A 7 "SAMTOOLS FLAGSTAT - DUPLICATES" ${summaryFile} | \
                tail -n 1)
    mapped=(${mapped//(/ }); mapped=${mapped[0]}; 
    duplicated=$(grep -A 6 "SAMTOOLS FLAGSTAT - DUPLICATES" ${summaryFile} | \
                tail -n 1)
    duplicated=(${duplicated/ / }); duplicated=${duplicated[0]}; 
    dupliPerce=$(echo "print((${duplicated}/${mapped})*100)"  | python)
    echo -e "\n% DUPLICATES" >> ${summaryFile}
    echo -e "${dupliPerce}\n" >> ${summaryFile}

    ## get % of reads aligning to spike-in
    # Only in case we look for a tag present in the spikes-in
    if grep -v "#" ${spikeInRef}.EpiListRenamed | grep -Fxq "${epiName}"; then
        # Get original chip name and look for counts to its spike-in
        # column 1 original name, column 2 in-house chip id
        epiOrigName=$(grep ${epiName} ${spikeInRef}.EpiListConversion | \
                        awk '{print $1}')
        lineNs=$(grep -n ${epiOrigName} ${spikeInRef}.fa | \
                cut --delimiter=":" --fields=1)

        # Get counts for focus tag
        lineNs=$(for li in ${lineNs}; do echo $(((${li}+1)/2)); done)
        spikeR1=$(for li in ${lineNs}; do 
                sed "${li}q;d" ${basePath}/spikeQC/counts/${filename}_R1.txt;
            done | awk '{sum+=$1;} END{print sum;}')
        spikeR2=$(for li in ${lineNs}; do 
                sed "${li}q;d" ${basePath}/spikeQC/counts/${filename}_R2.txt;
            done | awk '{sum+=$1;} END{print sum;}')
        spikeAdd=$((${spikeR1} + ${spikeR2}))
        echo -e "\nSpike focus counts" >> ${summaryFile}
        echo -e "${spikeAdd}\n" >> ${summaryFile}  

        # Get counts for all spike-ins
        spikeR1=$(cat ${basePath}/spikeQC/counts/${filename}_R1.txt | \
                awk '{sum+=$1;} END{print sum;}')
        spikeR2=$(cat ${basePath}/spikeQC/counts/${filename}_R2.txt | \
                awk '{sum+=$1;} END{print sum;}')
        spikeAddAll=$((${spikeR1} + ${spikeR2}))
        echo -e "Spike ALL counts" >> ${summaryFile}
        echo -e "${spikeAddAll}\n" >> ${summaryFile}  

    else
        echo -e "\nSpike focus counts" >> ${summaryFile}
        echo -e "TagNotInCUTANA\n" >> ${summaryFile}  

        # Get counts for all spike-ins
        spikeR1=$(cat ${basePath}/spikeQC/counts/${filename}_R1.txt | \
                awk '{sum+=$1;} END{print sum;}')
        spikeR2=$(cat ${basePath}/spikeQC/counts/${filename}_R2.txt | \
                awk '{sum+=$1;} END{print sum;}')
        spikeAddAll=$((${spikeR1} + ${spikeR2}))
        echo -e "Spike ALL counts" >> ${summaryFile}
        echo -e "${spikeAddAll}\n" >> ${summaryFile}  
    fi
    
    echo -e "mark duplicates - done ------------------------------------------------ \n"
    # store stage control info
    echo "DupCheck" >> ${stepControl}
else
    echo -e "mark duplicates - already done before ------------------------------------ \n"

fi 


echo -e "DATE CHECK \n $(date) \n" 



###################################
# qualiFlt: Qualimap + Filter by mapping q < 30
#################################
echo -e "Starting qualimap and Q<30 filtering -------------------------------\n"

bamSTR="${bamSTR}.q30"
bamSortMarkDupQ30="${EDITED_DIR}/BAM/${filename}.${bamSTR}.bam"

# check content of seventh line of step control file
linec=`sed "7q;d" ${stepControl}`
if [[ ${linec} != "qualiFlt" ]]; then 
    # Exception in thread "main" java.awt.AWTError: Can't connect to X11 
    #   window server...
    # I had to update qualimap script at line 41 to include 
    #   "-Djava.awt.headless=true" at java_options variable

    ${qualimap} bamqc --java-mem-size=${adjustedMem_Gb}G \
                -bam ${bamSortFlt} -gd MOUSE \
                -outdir ${basePath}/QC/qualimap/${filename}/ \
                -nt $SLURM_CPUS_PER_TASK \
                -outfile ${filename}.qualimap.report -outformat html

    echo -e "READS Q>30" >> ${summaryFile}
    samtools view -q 30 -c ${bamSortFlt} >> ${summaryFile}
    samtools view -q 30 -b ${bamSortFlt} > ${bamSortMarkDupQ30}

    if [[ $removeTemp == 'yes' ]] ; then
        rm ${bamSort}
    fi
    echo -e "Remove blacklist regions - done -----------------------------------------\n"
    # store stage control info
    echo "qualiFlt" >> ${stepControl}
else
    echo -e "Remove blacklist regions - already done before ------------------------------\n"
fi



# Remove chrM, chrUn, _random ... (remove sort.rmblackls.bam)
echo -e "Starting remove chrM and useless chromosomes ------------------------------\n"

bamSTR="${bamSTR}.rmchr"
bamSortMarkDupBlackChr="${EDITED_DIR}/BAM/${filename}.${bamSTR}.bam"

# check content of eigth line of step control file
linec=`sed "8q;d" ${stepControl}`
if [[ ${linec} != "Remove" ]]; then 
	samtools view -h ${bamSortMarkDupQ30} | \
    awk '(!index($3, "random")) && (!index($3, "chrUn")) && ($3 != "chrM") && ($3 != "chrEBV")' | \
    samtools view -Sb - > ${bamSortMarkDupBlackChr}
    samtools index ${bamSortMarkDupBlackChr} -@ $SLURM_CPUS_PER_TASK

    # QC: Show final reads
    echo -e "\nSAMTOOLS FLAGSTAT - FINAL READS" >> ${summaryFile}
    samtools flagstat ${bamSortMarkDupBlackChr} >> ${summaryFile}
    echo -e "\n" >> ${summaryFile}

    echo -e "Remove chrM and useless chromosomes - done ----------------------------------\n"
    # store stage control info
    echo "Remove" >> ${stepControl}
    if [[ $removeTemp == 'yes' ]] ; then
        rm ${bamSortMarkDupQ30}
    fi
else
    echo -e "Remove chrM and useless chromosomes - already done before--------------------------\n"
fi

# awk '(!index($1, "random")) && (!index($1, "chrUn")) && ($1 != "chrM") && ($1 != "chrEBV")' mm10.reordered.sizes


# Here i need to check if we get the number of mapped reads to assess if (ATAC numbers)
# > 50 million reads for open chrom and diff analysis
# > 200 million reads for TF footprinting based on empirical and computational stimations


#####################
# Trim Tn5 adaptors (+ strand 4bp to the right and - strand 5bp to the left)
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3959825
#https://training.galaxyproject.org/training-material/topics/epigenetics/tutorials/atac-seq/tutorial.html#check-insert-sizes
#####################

echo -e "Starting Tn5 adaptor removal ------------------------------\n"

if [ ! -e ${bamsPath} ]; then
    mkdir -p ${bamsPath}
fi

bamSTR="${bamSTR}.Tn5"
bamMarkDupBlackChrTn5="${EDITED_DIR}/BAM/${filename}.${bamSTR}.unsort.bam"
bamSortMarkDupBlackChrTn5="${bamsPath}/${filename}.${bamSTR}.bam"


# However, for peak calling, shifting of reads is not likely very important, as 
# it is a pretty minor adjustment and peaks are 100s of basepairs. The shifting 
# is only crucial when doing things where the exact position of the insertion 
# matters at single base resolution, e.g. TF motif footprinting.

# check content of ninth line of step control file
linec=`sed "9q;d" ${stepControl}`
if [[ ${linec} != "Tn5Adjust" ]]; then 
    # # use --ATACshift
    

    alignmentSieve --numberOfProcessors $SLURM_CPUS_PER_TASK \
            --ATACshift --bam ${bamSortMarkDupBlackChr} -o ${bamMarkDupBlackChrTn5}

    # the bam file needs to be sorted again
    samtools sort -@ $SLURM_CPUS_PER_TASK -o ${bamSortMarkDupBlackChrTn5} ${bamMarkDupBlackChrTn5}
    samtools index ${bamSortMarkDupBlackChrTn5} -@ $SLURM_CPUS_PER_TASK

    # samtools index -@ 8 sample1.shifted.bam
    rm ${bamMarkDupBlackChrTn5}
    echo -e "Tn5 adaptor removal - Finished -----------------------------\n"

    echo "Tn5Adjust" >> ${stepControl}

else
    echo -e "Tn5 adaptor removal  - already done before -----------------------------\n"
fi



###################################
# BigWigs normalised by spike in or cpm
#################################
# For spike-in normalization, a scale factor was calculated for each 
# sample by dividing the percent of total reads aligned to human genome 
# by the percent of total reads aligned to the spike-in barcodes (Scale 
# Factor = % Human Reads / % Spike-in Reads) and applying this factor 
# to adjust the total sequencing reads of each respective sample

# Normalization was done by taking the union set of all peaks from the 
# two H3K27ac ChIP samples (Hdac3-KO and Hdac3-WT), calculating the read 
# depth-normalized ratio of reads at each peak location, and then applying 
# the normalization factor derived from the spike-in panel barcodes to the 
# sequencing results - PMID: 32374402

# We first converted paired-end reads to fragments, then calculated the 
# enrichment of fragments mapping to each spike-in sequence after ChIP:
#     Esi=(barcode fragments in ChIP/barcode fragments in input).
# Next, we calculated H3K4me3 ChIP-Seq fragment coverage in all 25bp 
# nonoverlapping windows in the genome. The per-locus enrichment was 
# then calculated as:
#     Elocus=(fragment coverage in ChIP)/(fragment coverage in input).
# The histone modification density was then calculated as:
#     HMD(%)=100×(𝐸locus/𝐸H3K4me3_A).
# EH3K4me3_A is the Esi value for the H3K4me3_A spike_in. For our final 
# estimate of H3K4me3 signal, we normalized the HMD by the global average 
# HMD. The ChIP-Seq signal at peaks and other genomic intervals was 
# calculated as the sum of the signal in all overlapping 25bp windows.
# The K-MetStat panel cannot be used to normalize data from H3K9ac ChIP-Seq.
# From prev line i get that you can't normalise when you tag for a feature
# that is not contained in the spike-in panel? YOU CAN'T
# -PMID: 31444359

# deeptools bamCoverage
echo -e "Starting BigWigs --------------------------------------------------\n"

if [ ! -e ${basePath}/BigWig/ ]; then
	mkdir ${basePath}/BigWig/
fi

if [ ! -e ${EDITED_DIR}/BigWig/ ]; then
	mkdir ${EDITED_DIR}/BigWig/
fi


# Check if we can use the spike-in norm
spikeInScaled="no"
if grep -v "#" ${spikeInRef}.EpiListRenamed | grep -Fxq "${epiName}"; then
    spikeInScaled="yes"
    echo "Scaling Biwig by spike-in"
    bigWigOut="${basePath}/BigWig/${filename}.${bamSTR}.spike.bw"
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

    # get % of reads aligning to spike-in
    spikeAdd=$(grep -A 1 "Spike focus counts" ${summaryFile} | tail -n 1)
    spikePerce=$(echo "print((${spikeAdd}/${allRead})*100)"  | python)
else
    bigWigOut="${basePath}/BigWig/${filename}.${bamSTR}.CPM.bw"
fi


# check content of tenth line of step control file
linec=`sed "10q;d" ${stepControl}`
if [[ ${linec} != "BigWnorm1" ]]; then 

    if [[ ${spikeInScaled} == "yes" ]]; then
        # Factor = % Human Reads / % Spike-in Reads)
        #scale_factor=$(echo "${genomPerce} / $spikePerce" | bc -l)
        scale_factor=$(echo "print(${genomPerce}/${spikePerce})"  | python)
        echo -e "SCALING FACTOR\n$scale_factor\n" >> ${summaryFile}

        bamCoverage --binSize 5 --scaleFactor ${scale_factor}  \
            -b ${bamSortMarkDupBlackChrTn5} -of bigwig \
            -o ${basePath}/BigWig/${filename}.${bamSTR}.spike.bw \
            --numberOfProcessors $SLURM_CPUS_PER_TASK
    else
        echo -e "SCALING FACTOR\nCPM\n" >> ${summaryFile}
        
        
    fi
    bamCoverage --binSize 5 --normalizeUsing CPM --exactScaling \
            -b ${bamSortMarkDupBlackChrTn5} -of bigwig \
            -o ${basePath}/BigWig/${filename}.${bamSTR}.CPM.bw \
            --numberOfProcessors $SLURM_CPUS_PER_TASK

    
    # bedtools genomecov \
    #         -ibam $intervals \
    #         -bg -scale $scale_factor \
    #         > ${prefix}.${extension}

    # bedtools genomecov -bg -scale $scale_factor -i <sample>-converted.bed -g 
    # hg19.chrom.sizes > <sample>.bedGraph

    # bamCoverage \
    # --bam $input \
    # $args \
    # --scaleFactor ${scale} \
    # --numberOfProcessors ${task.cpus} \
    # --outFileName ${prefix}


	# bamCoverage --binSize 1 --normalizeUsing CPM --exactScaling \
    # -b ${bamSortMarkDupBlackChrTn5} -of bigwig \
    # -o ${bigWigOut} --numberOfProcessors $SLURM_CPUS_PER_TASK
	# bamCoverage --binSize 20 --normalizeUsing RPKM --effectiveGenomeSize 2913022398 \
	# -b ${BaseFolder}BAM/${filename}.sort.rmblackls.rmchr.bam -of bigwig \
	# -o ${BaseFolder}BigWig/${filename}.sort.rmblackls.rmchr.norm.bw
    
    echo -e "BigWig norm 1 - done ---------------------------------------------\n"
    # store stage control info
    echo "BigWnorm1" >> ${stepControl}
else
    echo -e "BigWig norm 1 - already done before ------------------------------\n"
fi



# I might remove this from the pipeline in the future
# echo -e "Starting BigWigs different normalization---------------------------\n"

# bigWigOut2="${EDITED_DIR}/BigWig/${filename}.sort.rmblackls.rmchr.gencovnorm.bw"

# # check content of eleventh line of step control file
# linec=`sed "12q;d" ${stepControl}`
# if [[ ${linec} != "BigWnorm2" ]]; then 
#     samtools view -b ${bamSortMarkDupBlackChrTn5} | \
#             genomeCoverageBed -ibam stdin -g $chr_genome_size -bg | \
#             $wigToBigWig -clip stdin $chr_genome_size ${bigWigOut2}
#     #samtools view -b ${bamSortMarkDupBlackChrTn5} | genomeCoverageBed -ibam stdin -g $chr_genome_size -bg > ${basePath}/BigWig/${filename}.sort.rmblackls.rmchr.gencovnorm.bedGraph
 
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

echo -e "PIPELINE END ---------------------------------------------------------"

seff $SLURM_JOBID

exit 0




# If control samples are provided in the sample sheet by default they will be 
# used to normalise the called peaks against non-specific background noise. 
# Control normalisation can be disabled using --use_control. Additionally it 
# may be necessary to scale control samples being used as background, 
# especially when read count normalisation methods have been used at earlier 
# stages in the pipeline.