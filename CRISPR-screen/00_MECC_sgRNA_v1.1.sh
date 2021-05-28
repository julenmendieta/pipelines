#!/bin/bash
# -*- ENCODING: UTF-8 -*-

#####===============================================================================
####   MECC_sgRNA_v1.sh
####   Counts for CrisprScreens
####   Maren Calleja @Cluster
####   Project: Crispr screens 
####   Wet-lab part: Ainhoa, David
####   Blah, blah
#####===============================================================================

##===============================================================================
## SLURM VARIABLES
#SBATCH --job-name=CRISPR
#SBATCH --cpus-per-task=6
#SBATCH --mem=10G
#SBATCH --time=05-10:00:00
#SBATCH -p short
#SBATCH --mail-type=END
#SBATCH --mail-user=mcallejac@unav.es
#SBATCH -o ./CRISPR_%a.log  # Standard output A = array ID, a = task ID
#SBATCH -e ./CRISPR_%a.err # Standard error
#find /home/mcallejac/David_Seq/data/V300042238_L02 -type f -name "WT*" -exec ln -s {} . ';'
#find /home/mcallejac/David_Seq/data/V300042238_L02 -type f -name "Cas9*" -exec ln -s {} . ';'
##sbatch --array=0-11%4 MECC_sgRNA_v1_495.sbs

# HOW TO RUN ME
# for i in *fastq.gz; do echo $i | sed 's/_R..fastq.gz//g' ; done | sort | uniq > samplesNames.txt
# N=`cat /home/jmendietaes/data/2021/CRISPR/sequencedData/merge4_492/samplesNames.txt | wc -l`
# sbatch --array=1-${N} 00_MECC_sgRNA_v1.1.sh \
#/home/jmendietaes/data/2021/CRISPR/sequencedData/merge4_492 \
#/home/jmendietaes/data/2021/CRISPR/allProcessed \
#/home/jmendietaes/referenceGenomes/mm10_reordered/mm10.reordered

##===============================================================================
## GLOBAL VARIABLES
# the base ones to modify

# -- base path were fastqs are located
PARENT_DIR="/home/jmendietaes/data/2021/CRISPR/sequencedData/merge4_492"
# -- regexp filter to extract the 20 nucleotide guides
# bulk library design (More info about pattern in line 54)
#flanquingSeq='CACCG(.{20})GTTTTAGAGC'
flanquingSeq='CACCG(.{20})GT{2,4}AGAGC'
# Alternative filter allowing from 2 to 4 T in the end 'CACCG(.{20})GT{2,4}AGAGC'
# new library (scRNAseq)
#flanquingSeq='CACCG(.{20})GTTTAAGAGC'




# these ones shouldnt need to be modified if we are using the right folder structure
FASTQ_DIR=$PARENT_DIR"/fastq"
EDITED_DIR=$PARENT_DIR"/pipelineOut"
#FASTQC_DIR=$EDITED_DIR"/fastQC"

## Guide REFERENCE 
## Comment by David These are the indexes for the alignment of the sgRNAs. But we will have just one genome index for every library.
# bowtie 1
GenomeIndex_all='/home/jmendietaes/referenceGenomes/sgRNA_indexes/bowtie1/finalGuides'
# bowtie 2
#GenomeIndex_all='/home/jmendietaes/referenceGenomes/sgRNA_indexes/bowtie2/finalGuides'
# change module and mapping code if you swich bowtie version in here

# path to the python script to extract the guide sequence
# the script will rely on the first '(.{' symbol in the pattern to get the number of nucleotides
#left-surrounding our set of nucleotides of interest
# the script will also assume that the number of nucleotides of interest is enclosed between '(.{'
#and '})', and get accordingly the (.{n}) nucleotides after the nucleotides preceding '(.{' in the matched string
# Hence, adding a filter pattern with '(.{' at the left of when we ask for our n nucleotides will cause fail
# first input is GZIP-ed fastq file path, then you can optionally add a regexp pattern 
#(default is 'CACCG(.{20})GTTTTAGAGC') and 'True' if you want the output fastq file
#to be also GZIP-ed
extractScript="/home/jmendietaes/programas/PhD/CRISPR-screen/00_NR_CRISPR-extract.py"

outGZ="False"
##===============================================================================
## Required Software
#module load FastQC/0.11.8-Java-1.8
#module load MultiQC/1.7-foss-2018b-Python-2.7.15
module load Bowtie/1.2.2-foss-2018b
#module load Bowtie2/2.3.4.2-foss-2018b
module load SAMtools/1.9-foss-2018b
#module load Perl/5.28.0-GCCcore-7.3.0
##===============================================================================

export PS4='$LINENO+ '
set -x

echo Maren Calleja  
echo me.callejac@gmail.com  mcallejac@unav.es
echo updated
echo julenmendieta92@gmail.com  jmendietaes@unav.es
echo Date ; date
echo script -- aln -- counts
echo LibraryType=SE
echo genome = mm10

# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

echo -e $SLURM_JOB_NODELIST 
#https://slurm.schedmd.com/faq.html#task_prolog
echo "print =========================================="
echo "print SLURM_JOB_ID = $SLURM_JOB_ID"
echo "print SLURM_JOB_NODELIST = $SLURM_JOB_NODELIST"
echo "print =========================================="

##===============================================================================
## Create or move around the corresponding directories
START_TIME=$SECONDS

if [ ! -e ${EDITED_DIR} ]; then
    mkdir -p ${EDITED_DIR}
fi

if [ ! -e ${EDITED_DIR}/QC/ ] && echo exists ; then
    mkdir -p ${EDITED_DIR}/QC/
fi

cd $EDITED_DIR
##===============================================================================
##Choose files# If we always get in list all files from folder
#FILES=($(ls $FASTQ_DIR | grep "fastq.gz" | sed -r 's/.{16}$//' | uniq))
# If we want to do only the files in samplesNames.txt file. Perfecto for re-running
#dead jobs
FILES=($(cat $FASTQ_DIR/samplesNames.txt))
filename=${FILES[$SLURM_ARRAY_TASK_ID]}

# get some paths
read1_path="${FASTQ_DIR}/${filename}_R1.fastq.gz"
stepControl="${EDITED_DIR}/QC/pipelineStep_${filename}.txt"
summaryFile="${EDITED_DIR}/QC/summary_${filename}.txt"

##===============================================================================
#############
## BEGIN: create summary file
#############

echo -e "Starting Summary file -------------------------------------- \n"

# check content of first line of step control file
linec=`sed "1q;d" ${stepControl}`
if [[ ${linec} != "Summary" ]]; then 
    echo -e "STARTING \n $(date) \n" 
    echo "SAMPLE: ${filename}" 
    echo -e "sample name\tfastq name\tread count\tmillions" >> ${summaryFile}

    # QC: read counts if file is gziped or not
    if file --mime-type ${read1_path} | grep -q gzip$; then
        Counts1="$(zcat ${read1_path} | echo $((`wc -l`/4)))"
    else
        echo "Not gzipped files"
        echo $read1_path
        #Counts1="$(cat ${RAW_FASTQ_DIR}/${filename}_R1_001.fastq | echo $((`wc -l`/4)))"
        #Counts2="$(cat ${RAW_FASTQ_DIR}/${filename}_R2_001.fastq| echo $((`wc -l`/4)))"
    fi

    echo -e "READ COUNTS \n" >> ${summaryFile}
    rc=$((Counts1/1000000))
    echo -e "${filename} \t ${filename}_R1 \t ${Counts1} \t ${rc}" >> ${summaryFile}
    
    echo -e "Summary file - done -------------------------------------- \n"
    # store stage control info
    echo "Summary" > ${stepControl}
else
    echo -e "Summary file - already done before -------------------------------------- \n"
fi


#############
# QC: fastqc
#############

# echo -e "Starting FastQC -------------------------------------- \n"

# if [ ! -e ${FASTQC_DIR} ]; then
#     mkdir -p ${FASTQC_DIR}
# fi

# # check content of second line of step control file
# linec=`sed "2q;d" ${stepControl}`
# if [[ ${linec} != "QC" ]]; then 
#     time fastqc ${read1_path} -o ${FASTQC_DIR} -t $SLURM_CPUS_PER_TASK

#     # MultiQC to join QC into a single HTML
#     # multiqc $fastq -o $fastq

#     echo -e "FastQC done -------------------------------------- \n"
#     # store stage control info
#     echo "QC" >> ${stepControl}
# else
#     echo -e "FastQC - already done before -------------------------------------- \n"
# fi


#############
# UZ: Unzip fastq files
#############
#echo -e "Unziping fasta files ------------------- \n"

#fastaUnz=${EDITED_DIR}/fastq_extr/${filename}_R1_001.fastq

# Create output dir
if [ ! -e ${EDITED_DIR}/fastq_extr/ ]; then
    mkdir -p ${EDITED_DIR}/fastq_extr/
fi

# check content of third line of step control file
#linec=`sed "3q;d" ${stepControl}`
#if [[ ${linec} != "UZ" ]]; then 
#    if file --mime-type ${read1_path} | grep -q gzip$; then
#        #gzip -d -c ${read1_path} > ${file/\.gz/};
#        gzip -d -c ${read1_path} > ${fastaUnz};
#    else 
#        echo "skipping $filename, not gzip"
#        exit 1
#    fi
#    echo -e "Unziping - done ------------------------ \n"

#    # store stage control info
#    echo "UZ" >> ${stepControl}
#else
#    echo -e "Unziping - already done before ------------------- \n"
#    
#fi

#############
# Extract: Extraction of guides sequence
#############

# Extraction
echo -e "Starting CRISPR-extract-carpools------------------- \n"

# Create output dir
if [[ ${outGZ} == "False" ]]; then
    extr_fastq="${EDITED_DIR}/fastq_extr/${filename}_extracted.fastq"
else
    extr_fastq="${EDITED_DIR}/fastq_extr/${filename}_extracted.fastq.gz"
fi

# Python script arguments
#		[Regular expression used to extract 20 bp target sequence]
#		[FastqFile] compressed, must finish in .gz
#		[True to GZIP reulting fastq]


# check content of second line of step control file
linec=`sed "2q;d" ${stepControl}`
if [[ ${linec} != "Extract" ]]; then 
    python ${extractScript} --fastq_path ${read1_path} --pattern ${flanquingSeq} --outgz ${outGZ}
    #perl ''${extractScript}'' ${flanquingSeq} ${fastaUnz} FALSE ${readId} >> ${EDITED_DIR}/QC/${filename}CRISPR-extract-carpools.log 2>&1

    echo -e "CRISPR-extract-carpools - done -------------------------------------- \n"
    # store stage control info
    echo "Extract" > ${stepControl}
else
    echo -e "CRISPR-extract-carpools - already done before ---------------------- \n"
fi

## David comment: Here, we use the flanking sequences to each sgRNA to extract the 20bp sgRNA sequence that will be used for the mapping process.
## Please notice that: In the bulk library design this is: CACCG(.{20})GTTTTAGAGC
## In the new library design (used for scRNAseq) this is: CACCG(.{20})GTTTAAGAGC


#############
# Alignemnt
#############

# Alignment and postprocess

echo -e "Starting Alignment -------------------------------------- \n"

## David comment: This is the alignment process with Bowtie2, usually we are above 90% aligment rate, it would be good to flag the samples that do not yield such %aligmmnet

if [ ! -e ${EDITED_DIR}/BAM/ ]; then
    mkdir -p ${EDITED_DIR}/BAM/
fi



# ESTO LO TENGO QUE CAMBIAR PARA QUE MAPEE EN UN INDEX CON TODOS LOS GUIDES Y LUEGO MEDIR CUANTOS GUIDES 
# SE HAN MAPEADO CON RESPECTO AL DE INTERES

samPath="${EDITED_DIR}/BAM/${filename}.sam"
# check content of third line of step control file
linec=`sed "3q;d" ${stepControl}`
if [[ ${linec} != "Align" ]]; then 
    # For relatively short reads (e.g. less than 50 bp) Bowtie 1 is sometimes faster and/or more sensitive.
    # http://bowtie-bio.sourceforge.net/bowtie2/faq.shtml
    bowtie2 -p $SLURM_CPUS_PER_TASK -x $GenomeIndex -U ${extr_fastq} -S ${samPath} >> ${EDITED_DIR}/BAM/${filename}bowtie2.log 2>&1
    
    # store stage control info
    echo "Align" >> ${stepControl}

    echo -e "Alignment - done -------------------------------------- \n"
else
    echo -e "Alignment - already done before -------------------------------------- \n"
 
fi


#############
# toBam: convert to BAM and create sorted index
#############

echo -e "Starting SAM to BAM -------------------------------------- \n"

bamPath="${EDITED_DIR}/BAM/${filename}.bam"
bamSortPath="${EDITED_DIR}/BAM/${filename}.sort.bam"
# check content of forth line of step control file
linec=`sed "4q;d" ${stepControl}`
if [[ ${linec} != "toBam" ]]; then 

    samtools view -o ${bamPath} -bhS -@ $SLURM_CPUS_PER_TASK ${samPath}
    # samtools view ${EDITED_DIR}/BAM/${filename}.bam| cut -f3 | sort | uniq -c > ${EDITED_DIR}/BAM/${filename}.cnt #After, I will discard this and retrieve indxstats ONLY
    # cat ${EDITED_DIR}/BAM/${filename}.cnt| awk '{print $1;}' > file2 # I think I will ddiscard this part and retrieve the .cnt only
    # cat ${EDITED_DIR}/BAM/${filename}.cnt| awk '{print $2;}' > file1
    # paste file1 file2 > ${EDITED_DIR}/BAM/${filename}.txt
    # rm file1 file2
    # head=( ${EDITED_DIR}/BAM/${filename}.txt COUNTS )
    # ( IFS=$'\t'; echo "${head[*]}"; cat ${EDITED_DIR}/BAM/${filename}.txt ) > ${EDITED_DIR}/BAM/${filename}.final.txt

    samtools sort -o ${bamSortPath} ${bamPath} 

    samtools index -b ${bamSortPath}

    # store stage control info
    echo "toBam" >> ${stepControl}

    echo -e "SAM to BAM - done -------------------------------------- \n"
else
    echo -e "SAM to BAM - already done before ---------------------------- \n"

fi


#############
# lastChecks: create flagstats and indxstats files and remove redundant data
#############

echo -e "Starting last checks -------------------------------------- \n"

# check content of fifth line of step control file
linec=`sed "5q;d" ${stepControl}`
if [[ ${linec} != "lastChecks" ]]; then 

    echo -e "SAMTOOLS FLAGSTAT \n" >> ${summaryFile}
    samtools flagstat ${bamSortPath} >> ${summaryFile}
    echo -e "\n"  >> ${summaryFile}

    # Indxstats = counts in this case
    samtools idxstats ${bamSortPath} > ${EDITED_DIR}/QC/${filename}.indxstats

    # Remove intermediate files, pickup the right ones... as needed
    rm -rf ${samPath}
    rm -rf ${bamPath}

    # store stage control info
    echo "lastChecks" >> ${stepControl}

    echo -e "last checks - done -------------------------------------- \n"
else
    echo -e "last checks - already done before ---------------------------- \n"

fi

echo -e "FINISHED... ------------------------------------------------------\n"

echo -e  Elapsed time $(($SECONDS - $START_TIME)) s

seff $SLURM_JOBID
