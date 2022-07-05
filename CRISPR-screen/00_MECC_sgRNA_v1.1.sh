#!/bin/bash
# -*- ENCODING: UTF-8 -*-

#####===============================================================================
####   MECC_sgRNA_v1.sh
####   Counts for CrisprScreens
####   Maren Calleja @Cluster
####   Updated: Julen Mendieta
####   Project: Crispr screens 
####   Wet-lab part: Ainhoa, David
####   Blah, blah
#####===============================================================================

##===============================================================================
## SLURM VARIABLES
#SBATCH --job-name=CRISPR
#SBATCH --cpus-per-task=6
#SBATCH --mem=10G
#SBATCH --time=00:30:00
#SBATCH -p short
#SBATCH -o /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.out  
#SBATCH -e /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.err 

##SBATCH --mail-type=END
##SBATCH --mail-user=mcallejac@unav.es
#find /home/mcallejac/David_Seq/data/V300042238_L02 -type f -name "WT*" -exec ln -s {} . ';'
#find /home/mcallejac/David_Seq/data/V300042238_L02 -type f -name "Cas9*" -exec ln -s {} . ';'
##sbatch --array=0-11%4 MECC_sgRNA_v1_495.sbs

# HOW TO RUN ME
# for i in *_R1_001.fastq.gz; do echo $i | sed 's/_R1_001.fastq.gz//g' ; done | sort | uniq > samplesNames.txt
# N=`cat samplesNames.txt | wc -l`
# sbatch --array=1-${N} /home/jmendietaes/programas/PhD/CRISPR-screen/00_MECC_sgRNA_v1.1.sh \
#/home/jmendietaes/data/2021/CRISPR/sequencedData/merge4-492 \
#/home/jmendietaes/data/2021/CRISPR/allProcessed/merge4-492 \
#/home/jmendietaes/referenceGenomes/sgRNA_indexes/bowtie2 \
#'CACCG(.{20})GT'

# the files stored in samplesNames.txt are the ones that will be analysed
# for lib in ../../Lib*; do cd $lib/demux_fastq; N=`cat samplesNames.txt | \
# wc -l`; path1=`realpath $lib` ; path2=`basename $path1`; \
# path2="/home/jmendietaes/data/2021/CRISPR/allProcessed/$path2"; \
# echo "sbatch --array=1-${N} 00_MECC_sgRNA_v1.1.sh ${path1} ${path2} \
# /home/jmendietaes/referenceGenomes/sgRNA_indexes/bowtie2 'CACCG(.{20})GT'"; done

# Files must have the S[0-9]_R[12]_001.fastq.gz structure, otherwise R script in
# next step will fail with rbind message

##===============================================================================
## GLOBAL VARIABLES
# the base ones to modify

# -- base path were fastqs are located
PARENT_DIR=$1
#PARENT_DIR="/home/jmendietaes/data/2021/CRISPR/sequencedData/merge4-492"
final_dir=$2
#final_dir="/home/jmendietaes/data/2021/CRISPR/allProcessed/merge4-492"

## Guide REFERENCE 
## Comment by David These are the indexes for the alignment of the sgRNAs. But we will have just one genome index for every library.
## In the end we go back to one index per library, so here we pass the main folder with
# each of the IDs named ID.fa as in the guide ID stated in every sample
GenomeIndex_all=$3
# bowtie 1
#GenomeIndex_all='/home/jmendietaes/referenceGenomes/sgRNA_indexes/bowtie1'
# bowtie 2
#GenomeIndex_all='/home/jmendietaes/referenceGenomes/sgRNA_indexes/bowtie2'
# change module and mapping code if you swich bowtie version in here
# -- regexp filter to extract the 20 nucleotide guides
# bulk library design (More info about pattern in line 54)
flanquingSeq=$4
#flanquingSeq='CACCG(.{20})GT{2,4}AGAGC'
#flanquingSeq='CACCG(.{20})GTTTTAGAGC'
# Alternative filter allowing from 2 to 4 T in the end 'CACCG(.{20})GT{2,4}AGAGC'
# new library (scRNAseq)
#flanquingSeq='CACCG(.{20})GTTTAAGAGC'
#flanquingSeq='CACCG(.{20})GT{2,4}'
# UTILIZAR ESTE DE AQUI
#flanquingSeq='CACCG(.{20})GT'
# Si se secuencias junto con ChIP (75bp)
# flanquingSeq='CACCG(.{18})'

# if the files have 001 or not, we add the string or not
extraSTR="_001"
#extraSTR=""


#######################################################################################
# these ones shouldnt need to be modified if we are using the right folder structure
FASTQ_DIR=$PARENT_DIR"/demux_fastq"
EDITED_DIR=$PARENT_DIR"/pipelineOut"
#FASTQC_DIR=$EDITED_DIR"/fastQC"

# path to the location of my git repo
gitPath="/home/jmendietaes/programas/PhD"

# path to the python script to extract the guide sequence
# the script will rely on the first '(.{' symbol in the pattern to get the number of nucleotides
#left-surrounding our set of nucleotides of interest
# the script will also assume that the number of nucleotides of interest is enclosed between '(.{'
#and '})', and get accordingly the (.{n}) nucleotides after the nucleotides preceding '(.{' in the matched string
# Hence, adding a filter pattern with '(.{' at the left of when we ask for our n nucleotides will cause fail
# first input is GZIP-ed fastq file path, then you can optionally add a regexp pattern 
#(default is 'CACCG(.{20})GTTTTAGAGC') and 'True' if you want the output fastq file
#to be also GZIP-ed
extractScript="${gitPath}/CRISPR-screen/00_NR_CRISPR-extract.py"

# path to R script for final guide analysis
RscriptP="${gitPath}/CRISPR-screen/01a_MECC_sgRNA_pre_analysis_v2.r"

# path to the python mapping report script
reportScript="${gitPath}/CRISPR-screen/01b_mergedMappingReport.py"
outGZ="False"

##===============================================================================
## Required Software
#module load FastQC/0.11.8-Java-1.8
#module load MultiQC/1.7-foss-2018b-Python-2.7.15
#module load Bowtie/1.2.2-foss-2018b
module load Bowtie2/2.3.4.2-foss-2018b
module load SAMtools/1.9-foss-2018b
#module load Perl/5.28.0-GCCcore-7.3.0

# set python paths
export PATH="/home/jmendietaes/programas/miniconda3/bin:$PATH"
export PYTHONPATH=/home/jmendietaes/programas/miniconda3/bin/python3.8
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

if [ ! -e ${final_dir}/QC/ ] && echo exists ; then
    mkdir -p ${final_dir}/QC/
fi

cd $EDITED_DIR
##===============================================================================
##Choose files# If we always get in list all files from folder
#FILES=($(ls $FASTQ_DIR | grep "fastq.gz" | sed -r 's/.{16}$//' | uniq))
# If we want to do only the files in samplesNames.txt file. Perfecto for re-running
#dead jobs
FILES=($(cat $FASTQ_DIR/samplesNames.txt))
filename=${FILES[$SLURM_ARRAY_TASK_ID - 1]}
sampleRunName=`basename ${final_dir}`

# get library ID to which map this file
mapLib=(${filename//_/ })
mapLib=${mapLib[2]}
GenomeIndex="$GenomeIndex_all/$mapLib/$mapLib.fa"
GenomeIndex_allGuide="$GenomeIndex_all/allGuides/finalGuides.fa"

# ensure that we have this index
if [[ ! -f ${GenomeIndex} ]]
then
    echo "Inferred index file does not exist"
    echo ${GenomeIndex}
    exit 1
fi

# get some paths
read1_path="${FASTQ_DIR}/${filename}_R1${extraSTR}.fastq.gz"
stepControl="${EDITED_DIR}/QC/pipelineStep_${filename}.txt"
summaryFile="${final_dir}/QC/summary_${filename}.txt"

if [ ! -e ${stepControl} ] ; then
    touch ${stepControl}
fi

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
    echo -e "READ COUNTS \n" >> ${summaryFile}
    echo -e "sample name\tfastq name\tread count\tmillions" >> ${summaryFile}

    # QC: read counts if file is gziped or not
    if file --mime-type ${read1_path} | grep -q gzip$; then
        Counts1="$(zcat ${read1_path} | echo $((`wc -l`/4)))"
    else
        echo "Not gzipped files"
        echo $read1_path
        exit 1
        #Counts1="$(cat ${RAW_FASTQ_DIR}/${filename}_R1_001.fastq | echo $((`wc -l`/4)))"
        #Counts2="$(cat ${RAW_FASTQ_DIR}/${filename}_R2_001.fastq| echo $((`wc -l`/4)))"
    fi

    rc=$((Counts1/1000000))
    echo -e "${sampleRunName} \t ${filename} \t ${Counts1} \t ${rc}" >> ${summaryFile}
    
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
    extr_fastq="${EDITED_DIR}/fastq_extr/${filename}_R1${extraSTR}_extracted.fastq"
    extr_unmap="${EDITED_DIR}/fastq_extr/${filename}_R1${extraSTR}_extracted_unmap.fastq"
else
    extr_fastq="${EDITED_DIR}/fastq_extr/${filename}_R1${extraSTR}_extracted.fastq.gz"
    extr_unmap="${EDITED_DIR}/fastq_extr/${filename}_R1${extraSTR}_extracted_unmap.fastq.gz"
fi

# Python script arguments
#		[Regular expression used to extract 20 bp target sequence]
#		[FastqFile] compressed, must finish in .gz
#		[True to GZIP reulting fastq]


# check content of second line of step control file
linec=`sed "2q;d" ${stepControl}`
if [[ ${linec} != "Extract" ]]; then 
    echo -e "\nGUIDE EXTRACTION\n" >> ${summaryFile}
    python ${extractScript} --fastq_path ${read1_path} --pattern ${flanquingSeq} --outgz ${outGZ} >> ${summaryFile}
    mv ${FASTQ_DIR}/${filename}_R1${extraSTR}_extracted.fastq ${EDITED_DIR}/fastq_extr/
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
samPath_unM="${EDITED_DIR}/BAM/${filename}.unMap.sam"
# check content of third line of step control file
linec=`sed "3q;d" ${stepControl}`
if [[ ${linec} != "Align" ]]; then 
    # For relatively short reads (e.g. less than 50 bp) Bowtie 1 is sometimes faster and/or more sensitive.
    # http://bowtie-bio.sourceforge.net/bowtie2/faq.shtml
    #bowtie -p $SLURM_CPUS_PER_TASK ${GenomeIndex} ${extr_fastq} -S ${samPath} --best

    # The mapping takes 22 secods, mapping percentaje difference goes from 96.58 in bowtie 2 to 97.23 in bowtie
    # but I like that bowtie 2 states number of alignments of reads in more than one ref
    echo -e "\nALIGNMENT\n" >> ${summaryFile}
    bowtie2 -p $SLURM_CPUS_PER_TASK -x $GenomeIndex -U ${extr_fastq} \
            -S ${samPath} --un ${extr_unmap} >> ${summaryFile} 2>&1
    
    # with --un we write unmapped reads to another fastq
    # then we try to map them to all the reference guides
    echo -e "\nALIGNMENT of previously unmapped\n" >> ${summaryFile}
    bowtie2 -p $SLURM_CPUS_PER_TASK -x $GenomeIndex_allGuide -U ${extr_unmap} \
            -S ${samPath_unM} >> ${summaryFile} 2>&1

    # Default mode: search for multiple alignments, report the best one
    # Bowtie 2 does not guarantee that the alignment reported is the best possible in terms of alignment score.
    # Increasing -R makes Bowtie 2 slower, but increases the likelihood that it will report the correct alignment 
    # for a read that aligns many places.
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

#samPathUniq="${EDITED_DIR}/BAM/${filename}.uniqM.sam"
bamPath="${EDITED_DIR}/BAM/${filename}.bam"
bamSortPath="${EDITED_DIR}/BAM/${filename}.sort.bam"

bamPath_unM="${EDITED_DIR}/BAM/${filename}.unMap.bam"
bamSortPath_unM="${EDITED_DIR}/BAM/${filename}.unMap.sort.bam"
# check content of forth line of step control file
linec=`sed "4q;d" ${stepControl}`
if [[ ${linec} != "toBam" ]]; then 
    # first we filter out reads mapping more than twice
    #samtools view -h ${samPath} | grep -v XS:i > ${samPathUniq}
    #skipped=`samtools view ${samPath} | grep XS:i | wc -l`
    skipped=0
    echo -e "\nFILTER READS MAPPING MORE THAN ONCE\n" >> ${summaryFile}
    echo -e "$skipped reads filtered out\n" >> ${summaryFile}
    echo -e "Number of reads that map more than once, and their top target\n" >> ${summaryFile}
    samtools view -h ${samPath} | grep XS:i | awk '{print $3}' | sort | uniq -c >> ${summaryFile}
    # Then we convert to bam
    samtools view -o ${bamPath} -bhS -@ $SLURM_CPUS_PER_TASK ${samPath}
    samtools view -o ${bamPath_unM} -bhS -@ $SLURM_CPUS_PER_TASK ${samPath_unM}
    # samtools view ${EDITED_DIR}/BAM/${filename}.bam| cut -f3 | sort | uniq -c > ${EDITED_DIR}/BAM/${filename}.cnt #After, I will discard this and retrieve idxstats ONLY
    # cat ${EDITED_DIR}/BAM/${filename}.cnt| awk '{print $1;}' > file2 # I think I will ddiscard this part and retrieve the .cnt only
    # cat ${EDITED_DIR}/BAM/${filename}.cnt| awk '{print $2;}' > file1
    # paste file1 file2 > ${EDITED_DIR}/BAM/${filename}.txt
    # rm file1 file2
    # head=( ${EDITED_DIR}/BAM/${filename}.txt COUNTS )
    # ( IFS=$'\t'; echo "${head[*]}"; cat ${EDITED_DIR}/BAM/${filename}.txt ) > ${EDITED_DIR}/BAM/${filename}.final.txt

    samtools sort -o ${bamSortPath} ${bamPath} 
    samtools index -b ${bamSortPath}

    samtools sort -o ${bamSortPath_unM} ${bamPath_unM} 
    samtools index -b ${bamSortPath_unM}

    # store stage control info
    echo "toBam" >> ${stepControl}

    echo -e "SAM to BAM - done -------------------------------------- \n"
else
    echo -e "SAM to BAM - already done before ---------------------------- \n"

fi


#############
# lastChecks: create flagstats and idxstats files and remove redundant data
#############

if [ ! -e ${final_dir}/idxstats/ ]; then
    mkdir -p ${final_dir}/idxstats/
fi

echo -e "Starting last checks -------------------------------------- \n"

# check content of fifth line of step control file
linec=`sed "5q;d" ${stepControl}`
if [[ ${linec} != "lastChecks" ]]; then 

    echo -e "\nSAMTOOLS FLAGSTAT\n" >> ${summaryFile}
    samtools flagstat ${bamSortPath} | head -5 | tail -1 >> ${summaryFile}
    echo -e "\nSAMTOOLS FLAGSTAT OF UNMAPPED\n" >> ${summaryFile}
    samtools flagstat ${bamSortPath_unM} | head -5 | tail -1 >> ${summaryFile}
    echo -e "\n"  >> ${summaryFile}

    # idxstats = counts in this case
    samtools idxstats ${bamSortPath} > ${final_dir}/idxstats/${filename}.idxstats
    samtools idxstats ${bamSortPath_unM} > ${final_dir}/idxstats/${filename}.unMap.idxstats

    # Remove intermediate files, pickup the right ones... as needed
    rm -rf ${samPath}
    rm -rf ${samPath_unM}
    rm -rf ${bamPath}
    rm -rf ${bamPath_unM}
    rm -rf ${extr_fastq}

    # store stage control info
    echo "lastChecks" >> ${stepControl}

    echo -e "last checks - done -------------------------------------- \n"
else
    echo -e "last checks - already done before ---------------------------- \n"

fi

#############
# Rscript: prepare code to run R script
#############

if [ ! -e ${final_dir}/RSession/ ]; then
    mkdir -p ${final_dir}/RSession/
fi

echo -e "Starting Rscript text ---------------------------------- \n"

# check content of sixth line of step control file
linec=`sed "6q;d" ${stepControl}`
if [[ ${linec} != "Rscript" ]]; then 
    echo "module load R/4.0.5-foss-2020b" > ${final_dir}/RSession/toRunR.txt
    echo "python ${reportScript} -ip ${final_dir}" >> ${final_dir}/RSession/toRunR.txt
    echo "Rscript --vanilla ${RscriptP} ${final_dir} > ${final_dir}/RSession/${sampleRunName}.Rout.txt" >> ${final_dir}/RSession/toRunR.txt
    echo "cd ${final_dir}/RSession" >> ${final_dir}/RSession/toRunR.txt
    echo "zip ${sampleRunName}.zip *" >> ${final_dir}/RSession/toRunR.txt

    # store stage control info
    echo "Rscript" >> ${stepControl}

    echo -e "Rscript - done -------------------------------------- \n"
else
    echo -e "Rscript - already done before ---------------------------- \n"
fi


echo -e "FINISHED... ------------------------------------------------------\n"

echo -e  Elapsed time $(($SECONDS - $START_TIME)) s

seff $SLURM_JOBID
