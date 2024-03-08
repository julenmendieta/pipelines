#!/bin/bash
# -*- ENCODING: UTF-8 -*-

##===============================================================================
## SLURM VARIABLES
#SBATCH --job-name=EMseq_fqToBw_Biscuit
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
# sbatch --array=1-${N} /home/jmendietaes/programas/pipelines/EM-seq/cluster/01_EMseq_fqTObw_allBiscuit.sh \
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
# https://huishenlab.github.io/biscuit/
# Remov efollowing if not used
# https://github.com/semenko/serpent-methylation-pipeline/tree/main
# https://github.com/nebiolabs/EM-seq/blob/master/em-seq.nf
# https://nf-co.re/methylseq/2.4.0/parameters#clip_r1


##===============================================================================
# Path to extra scripts
subScripts="/home/jmendietaes/programas/pipelines/EM-seq/cluster/sub-scripts"

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
module load SAMtools/1.12-GCC-10.2.0
module load BEDTools/2.27.1-foss-2018b
module load zlib/1.2.11-GCCcore-11.2.0
export PATH="/home/jmendietaes/programas/miniconda3/envs/DNAme/bin:$PATH"

bgzip=~/programas/miniconda3/bin/bgzip

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
#Aling # Mapping
############
echo -e "Starting Alignment -------------------------------------- \n"

if [ ! -e ${EDITED_DIR}/BAM/ ]; then
    mkdir -p ${EDITED_DIR}/BAM/
fi

bamFile="${EDITED_DIR}/BAM/${filename}.sort.markdup.bam"
bamFile_flt="${EDITED_DIR}/BAM/${filename}.sort.markdup.flt.bam"
bamFile_fltV="${EDITED_DIR}/BAM/${filename}.sort.markdup.fltV.bam"
# check content of third line of step control file
linec=`sed "2q;d" ${stepControl}`
if [[ ${linec} != "Align" ]]; then 

    # Align sequencing reads to the reference
    # Gzipped FASTQ files can also be used
    biscuit align -@ ${nCPU} ${REFERENCE_DIR}.fa ${read1_path} ${read2_path} | \
        dupsifter ${REFERENCE_DIR}.fa \
                -O ${EDITED_DIR}/QC/${filename}.output.dupsifter.stat | \
        samtools sort -@ ${nCPU} -o ${bamFile} -O BAM -
    samtools index ${bamFile}

    # Store duplicate stats
    dupsifterOut=${EDITED_DIR}/QC/${filename}.output.dupsifter.stat
    echo -e "\nDupsifter" >> ${summaryFile}
    grep "individual reads" ${dupsifterOut} >> ${summaryFile}
    grep "reads with both reads mapped" ${dupsifterOut}  >> ${summaryFile}
    grep "both reads marked as duplicate" ${dupsifterOut}  >> ${summaryFile}
    indivR=$(grep "individual reads" ${dupsifterOut}  | awk '{print $7}')
    alignedR=$(grep "reads with both reads mapped" ${dupsifterOut}  | awk '{print $9}')
    dupliR=$(grep "both reads marked as duplicate" ${dupsifterOut}  | awk '{print $11}')
    perceAlign=$(echo "print(round((${alignedR}/(${indivR}/2))*100,2))" | python)
    perceDupli=$(echo "print(round((${dupliR}/${alignedR})*100,2))" | python)
    echo "% of aligned reads: ${perceAlign}" >> ${summaryFile}
    echo "% of duplicated reads (from aligned): ${perceDupli}" >> ${summaryFile}

    # Optional
    # Filter Reads by Bisulfite Conversion
    # The -p flag outputs the counts in a table, instead of as a tag in the BAM file
    biscuit bsconv ${REFERENCE_DIR}.fa ${bamFile}  ${bamFile_flt}
    # Validate Bisulfite Conversion Label
    biscuit bsstrand ${REFERENCE_DIR}.fa ${bamFile_flt} ${bamFile_fltV}
    # f: OT/CTOT (BSW) strand
    # r: OB/CTOB (BSC) strand
    # c: conflicting strand information
    # u: unintelligible strand source (unknown)

        
    echo -e "Alignment - done -------------------------------------- \n"
    # store stage control info
    echo "Align" >> ${stepControl}

else
    echo -e "Alignment - already done before ------------------------------ \n"
fi


##############################
# Remove chrM, chrUn, _random ... (remove sort.rmdup.rmblackls.bam)
##############################

echo -e "Starting remove chrM and useless chromosomes ------------------------------\n"

if [ ! -e ${bamsPath} ]; then
    mkdir -p ${bamsPath}/allBiscuit
fi

bamFile_fltVChr="${bamsPath}/allBiscuit/${filename}.sort.markdup.fltV.rmchr.bam" 

# check content of eigth line of step control file
linec=`sed "3q;d" ${stepControl}`
if [[ ${linec} != "Remove" ]]; then 
	samtools view -h ${bamFile_fltV} | \
    awk '(!index($3, "random")) && (!index($3, "chrUn")) && ($3 != "chrM") && ($3 != "chrEBV")' | \
    samtools view -Sb - > ${bamFile_fltVChr}
    samtools index ${bamFile_fltVChr} -@ ${nCPU}

    # QC: Show final reads
    echo -e "\nSAMTOOLS FLAGSTAT - FINAL READS" >> ${summaryFile}
    samtools flagstat ${bamFile_fltVChr} >> ${summaryFile}
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
# Biscuit pileup
#################################

echo -e "Starting Biscuit pileup ------------------------------\n"

pileupP=${basePath}/biscuit/pileup/
if [ ! -e ${pileupP} ]; then
    mkdir -p "${pileupP}"
fi

# check content of eigth line of step control file
linec=`sed "4q;d" ${stepControl}`
if [[ ${linec} != "convrtd" ]]; then 
    
    # Create a pileup VCF of DNA methylation and genetic information
    # Also compresses and indexes the VCF
    # allows the user to compute cytosine retention and callable SNP mutations.
    # BISCUIT has the ability to put mutation calls and DNA methylation 
    # measurements from multiple samples next to each other in the output VCF 
    # by providing biscuit pileup with more than one input BAM
    biscuit pileup -@ ${nCPU} -o ${pileupP}/${filename}_my_pileup.vcf \
                    ${REFERENCE_DIR}.fa ${bamFile_fltVChr}
    ${bgzip} -@ ${nCPU} ${pileupP}/${filename}_my_pileup.vcf
    tabix -p vcf ${pileupP}/${filename}_my_pileup.vcf.gz

    # Extract DNA methylation into BED format
    # Also compresses and indexes the BED
    biscuit vcf2bed -t cg ${pileupP}/${filename}_my_pileup.vcf.gz \
                    > ${pileupP}/${filename}_methylation_data.bed
    ${bgzip} -@ ${nCPU} ${pileupP}/${filename}_methylation_data.bed
    tabix -p bed ${pileupP}/${filename}_methylation_data.bed.gz

    #  Merge neighboring C and G in CpG context
    biscuit mergecg ${REFERENCE_DIR}.fa \
                    ${pileupP}/${filename}_methylation_data.bed.gz \
                    1> ${pileupP}/${filename}_mergecg_data.bed
    ${bgzip} -@ ${nCPU} ${pileupP}/${filename}_mergecg_data.bed
    tabix -p bed ${pileupP}/${filename}_mergecg_data.bed.gz

    # Get file for Bismar COV (Methylation coverage track)
    biscuit vcf2bed -c ${pileupP}/${filename}_my_pileup.vcf.gz \
                    > ${pileupP}/${filename}_beta_m_u.bed
    # Chromosome
    # Start position (0-based)
    # End position
    # Methylation percentage
    # M (number of methylated reads covering locus)
    # U (number of unmethylated reads covering locus)
    awk -v OFS='\t' '{ print $1, $2+1, $3, $4, $5, $6 }' \
            ${pileupP}/${filename}_beta_m_u.bed \
            > ${pileupP}/${filename}_beta_m_u.cov
    # This file will keep methylation across the strands separate, as is normally 
    # done in biscuit vcf2bed. To merge methylation across strands, run:
    #biscuit vcf2bed my_pileup.vcf.gz | \
    #biscuit mergecg -c /path/to/my_reference.fa - | \
    #awk -v OFS='\t' '{ print $1, $2+1, $3-1, $4, $5, $6 }' > my_merged_beta_m_u.cov



    echo -e "mark of unconverted reads - done ----------------------------------\n"
    # store stage control info
    echo "convrtd" >> ${stepControl}

else
    echo -e "mark of unconverted reads - already done before--------------------------\n"
fi

#samtools view -F 512 ${samSortRmDupBlackChrConv}
#samtools view -o ${bamSortRmDup} -@ $SLURM_CPUS_PER_TASK -bh -F 512 ${bamSortMarkDup}

echo -e "END --------------------------------------------------"

seff $SLURM_JOBID

exit 0



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