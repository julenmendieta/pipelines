#!/bin/bash
# -*- ENCODING: UTF-8 -*-

##===============================================================================
## SLURM VARIABLES
#SBATCH --job-name=mRNA_countA
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --time=00-24:00:00
#SBATCH -p short
#SBATCH -o /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.out  
#SBATCH -e /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.err 
##SBATCH --dependency=afterany:793158


##SBATCH --mail-type=END
##SBATCH --mail-user=user@mail.es
# HOW TO RUN ME
# sbatch /home/jmendietaes/programas/pipelines/mRNA/cluster/01_RNAseq_countAnalysis.sh \
#/home/jmendietaes/data/2021/mRNA/allProcessed \
#/home/jmendietaes/referenceGenomes/mm10_reordered/STAR/mm10.reordered \
#01_projectRestart

# NOTES
# This pipeline is mainly based on:
#   https://github.com/nf-core/rnaseq
#   With few minor changes
# Also looked at:
#   https://bookdown.org/jean_souza/PreProcSEQ
#   http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#input-data

##==============================================================================
## GLOBAL VARIABLES

# path where we have the folder structure for chip analysis 
# (inside we have ${basePath}/bamfiles/valid/)
basePath=$1
#basePath=/home/jmendietaes/data/2021/mRNA/allProcessed

##  REFERENCE Genome info (we use it to catch bed file of same GTF as 00_...)
REFERENCE_DIR=$2
#REFERENCE_DIR="/home/jmendietaes/referenceGenomes/mm10_reordered/STAR/mm10.reordered"
indexOutP=$(dirname ${REFERENCE_DIR})
#geneBed=$(realpath ${indexOutP}/genes/*bed)
geneGTF=$(realpath ${indexOutP}/genes/*gtf)
# use singleCell GTF for test
#geneGTF=/home/jmendietaes/data/2021/singleCell/additionalFiles/refdata-gex-mm10-2020-A/genes/genes_mainGenes.gtf


## Project name
# We will find salmon output and STAR BAM files in there using our defined
#   folder structure
projectName=$3
#projectName=01_projectRestart

# Where to look for Salmon count files
# A path with folder named by sample and containing quant.sf inside
# I link the interest folders from allProcessed/counts/salmon to
#       another location to divide data by projects
#salmonC="${basePath}/countGroups/${projectName}"
# Were to look for bam files
bamsPath="${basePath}/bamfiles/wholeGenome/valid/${projectName}"
outpath=${basePath}"/furtherAnalysis/${projectName}"
salmonQin=${basePath}/counts/salmon/${projectName}

# Path to nextflow scripts
subScripts="/home/jmendietaes/programas/pipelines/mRNA/cluster/sub-scripts"


##==============================================================================
## Required Software

# Path to R
R="/home/jmendietaes/programas/miniconda3/envs/Renv/bin/Rscript"
# Path to Qualimap
qualimap="/home/jmendietaes/programas/pipelines/qualimap/qualimap_v2.2.1/qualimap"
# Java needed for Qualimap
module load Java/1.8.0_192
# Path to numfmt
numfmt=~/programas/miniconda3/bin/numfmt
# Variables from the job
nCPU=$SLURM_CPUS_PER_TASK

allbams=$(find ${bamsPath}/*bam -printf "${bamsPath}/%f\n" | \
        tr '\n' ' ')
allLabels=`for i in $allbams; do basename ${i} | cut -d '.' -f 1; done | \
        tr '\n' ' '`

##==============================================================================

# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

##==============================================================================


if [ ! -e ${outpath}/counts ]; then
    mkdir -p ${outpath}/counts
    mkdir -p ${outpath}/QC/biotypes
fi


if [ ! -e ${outpath}/QC/deseq2 ]; then
    mkdir -p ${outpath}/QC/deseq2
    mkdir -p ${outpath}/QC/biotypes
    mkdir -p ${outpath}/QC/qualimap
    mkdir -p ${outpath}/QC/dupradar
fi

##==============================================================================
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
        # If $1 is older than any $2
        if [[ $1 -ot ${tfile} ]] ; then
        analyse="yes"
        echo $1" older than "${tfile}
        fi
    done
    fi
}

# Load assigned RAM
adjustedMem=$(echo "print(int(round($SLURM_MEM_PER_NODE*0.96,0)/1000))" | python3)
# fast check in cases where RAM given in Gb (value will be lower than 1 in most cases)
# 100Gb / 1000 = 0.1
if [ "$adjustedMem" -lt "5" ]; then
    adjustedMem=$(echo "print(int(round($SLURM_MEM_PER_NODE*0.96,0)))" | python3)
fi  
echo  "Memory in Gb:"
echo ${adjustedMem}
# Convert to bytes
adjustedMem_bytes=$(printf ${adjustedMem} | ${numfmt} --from-unit=Mi)
adjustedMem_Gb=${adjustedMem}
echo  "Memory in bytes:"
echo ${adjustedMem_bytes}

##==============================================================================

##############################
# Transcript level abundance
##############################

# As Salmon recommends, aggregate transcript-level salmon abundance estimates 
#   to the gene level using tximport 
#   (https://bioconductor.org/packages/release/bioc/html/tximport.html)
# It is more versatile, exposes more features, and allows considering 
#   multi-sample information during aggregation
# Outputs:
#   -tpm: Scaled to library size
#   -length_scaled: scaled using the average transcript length over samples 
#       and then the library size 

salmonOut=${outpath}/counts
cd ${salmonOut}

# Check if we already have obtained gene counts
fileNotExistOrOlder "${salmonOut}/tximportMerge.transcript_counts.rds" \
            "${salmonQin}/*/quant.sf"
if [[ ${analyse} == "yes" ]]; then
    echo "Aggregate transcript-level salmon abundance estimates"
    # --salmon should be a folder with files names as samples and 
    #       containing Salmons quantification output quant.sf
    python ${subScripts}/salmon_tx2gene.py \
        --gtf ${geneGTF} \
        --salmon ${salmonQin} \
        --id "gene_id" \
        --extra "gene_name" \
        -o ${salmonOut}/salmon_tx2gene.tsv

    # Following script requires tximeta, but not sure for what
    ${R} ${subScripts}/salmon_tximport.r \
        Infer \
        ${salmonQin} \
        ${salmonOut}/tximportMerge

    # Store R rds sessions with summarized experiment of:
    #   Gene counts
    #   Gene counts scaled by gene length
    #   Gene counts scaled
    #   Transcript counts
    for f2 in gene_counts.tsv gene_counts_length_scaled.tsv \
    gene_counts_scaled.tsv; do
        ${R} ${subScripts}/salmon_summarizedexperiment.r \
            Infer \
            ${salmonOut}/tximportMerge.${f2} \
            ${salmonOut}/tximportMerge.gene_tpm.tsv
    done
    ${R} ${subScripts}/salmon_summarizedexperiment.r \
        Infer \
        ${salmonOut}/tximportMerge.transcript_counts.tsv \
        ${salmonOut}/tximportMerge.transcript_tpm.tsv
fi

# Get samples Biotype, and make sure is the same
samples=$(head -n 1 ${salmonOut}/tximportMerge.gene_tpm.tsv | \
                tr '\t' '\n' | \
                grep -v gene_ | sed 's/\./-/g')

lTypes=$(for s in ${samples}; do 
    grep -A 1 "LIBRARY TYPE" ${basePath}/QC/allStepsSummary/summary_${s}*txt |\
    tail -n 1; done | sort | uniq)
nTypes=$(echo ${lTypes} | wc -w)
if [ "$nTypes" -gt 1 ] ; then
    echo "ERROR: More than one library type"
    echo $nTypes 
    exit 1
fi


##############################
# QC: With deseq2
##############################

# Check if we have already have done the Deseq2 QC
fileNotExistOrOlder "${outpath}/QC/deseq2/deseq2.plots.pdf" \
            "${salmonOut}/tximportMerge.gene_counts_length_scaled.tsv"
if [[ ${analyse} == "yes" ]]; then
    ${R} ${subScripts}/deseq2_qc.r \
        --count_file ${salmonOut}/tximportMerge.gene_counts_length_scaled.tsv \
        --outdir ${outpath}/QC/deseq2/ \
        --cores ${nCPU} 
fi


###################################
# QC: Measure gene biotypes
#################################
# In a high quality data set, we expect the different biotypes to be equally
#   present. Deviations from this pattern may correspond to different 
#   biological groups in the data set, especially in case of different cell 
#   lines, but could also indicate batch effects or low quality samples.

header="# id: 'biotype_counts'\n\
# section_name: 'Biotype Counts'\n\
# description: 'shows reads overlapping genomic features of different biotypes,\n\
#     \tcounted by <a href='http://bioinf.wehi.edu.au/featureCounts'>featureCounts</a>.'\n\
# plot_type: 'bargraph'\n\
# anchor: 'featurecounts_biotype'\n\
# pconfig:\n\
#     \tid: 'featurecounts_biotype_plot'\n\
#     \ttitle: 'featureCounts: Biotypes'\n\
#     \txlab: '# Reads'\n\
#     \tcpswitch_counts_label: 'Number of Reads'"

#SUBREAD_FEATURECOUNTS 
# Define strandness
# 0 (unstranded), 1 (stranded) and 2 (reversely stranded)
lType_=0
if [[ ${lTypes} == *"SF" ]] ; then
    lType_=1
elif [[ ${lTypes} == *"SR" ]] ; then
    lType_=2
fi

for bam in ${allbams}; do
    sampleId=$(basename ${bam} | cut -d '.' -f 1)
    # Check if we have already have done the biotype QC with this BAM
    fileNotExistOrOlder "${outpath}/QC/biotypes/${sampleId}_biotype_counts_mqc.tsv" \
            "${bam}"
    if [[ ${analyse} == "yes" ]]; then

        featureCounts -p  \
            -T ${nCPU} \
            -g gene_type \
            -a ${geneGTF} \
            -s ${lType_} \
            -o ${outpath}/QC/biotypes/${sampleId}.featureCounts.txt \
            ${bam}


        echo -e $header > ${outpath}/QC/biotypes/${sampleId}_biotype_counts_mqc.tsv
        cut -f 1,7 ${outpath}/QC/biotypes/${sampleId}.featureCounts.txt | \
            tail -n +3 >> ${outpath}/QC/biotypes/${sampleId}_biotype_counts_mqc.tsv

        allBiotypes=$(grep -v "^#" \
            ${outpath}/QC/biotypes/${sampleId}_biotype_counts_mqc.tsv | \
            cut -f 1 | tr '\n' ' ')
        python ${subScripts}/mqc_features_stat.py \
            ${outpath}/QC/biotypes/${sampleId}_biotype_counts_mqc.tsv \
            -s ${sampleId} \
            -f ${allBiotypes} \
            -o ${outpath}/QC/biotypes/${sampleId}_biotype_counts_mqc.tsv
    fi
done
# Delete unnecesary files
rm ${outpath}/QC/biotypes/*featureCounts*
# Get biotype plot
python ${subScripts}/biotypePlot.py ${outpath}/QC/biotypes/
# Delete remaining files
#rm ${outpath}/QC/biotypes/*counts_mqc.tsv

# I can use DESeq scores to make these bigwigs
# Generate bigwig with --scaleFactor 1/<DESeq2's sizeFactor>


###################################
# QC: Qualimap
#################################
# http://qualimap.conesalab.org/doc_html/command_line.html
# Define strandness
# strand-specific-forward, strand-specific-reverse or non-strand-specific
lType_="non-strand-specific"
if [[ ${lTypes} == *"SF" ]] ; then
    lType_="strand-specific-forward"
elif [[ ${lTypes} == *"SR" ]] ; then
    lType_="strand-specific-reverse"
fi


# This step is quite slow
echo "Running Qualimap"
cd ${outpath}/QC/qualimap/
for bam in ${allbams}; do
    sampleId=$(basename ${bam} | cut -d '.' -f 1)

    # Check if we have already have run Qualimap with this BAM
    fileNotExistOrOlder "${outpath}/QC/qualimap/${sampleId}/qualimapReport.html" \
            "${bam}"
    if [[ ${analyse} == "yes" ]]; then
        echo ${sampleId}
        mkdir -p ${outpath}/QC/qualimap/${sampleId}
        ${qualimap} rnaseq --java-mem-size=${adjustedMem_Gb}G \
                -bam ${bam} \
                -gtf ${geneGTF} \
                -p ${lType_} \
                --paired \
                -outdir ${outpath}/QC/qualimap/${sampleId}
    fi
done
#QUALIMAP_RNASEQ


###################################
# QC: dupRadar
#################################
#DUPRADAR
# 0(unstranded)/1(forward)/2(reverse)
lType_=0
if [[ ${lTypes} == *"SF" ]] ; then
    lType_=1
elif [[ ${lTypes} == *"SR" ]] ; then
    lType_=2
fi

for bam in ${allbams}; do
    ${R} ${subScripts}/dupradar.r ${bam} \
        ${outpath}/QC/dupradar \
        ${geneGTF} \
        ${lType_} \
        paired \
        ${nCPU}
done


###################################
# Diff: Differential expression analysis
#################################
# The DESeq2 model internally corrects for library size, so transformed or 
#       normalized values such as counts scaled by library size should not 
#       be used as input.
# Note that the tximport-to-DESeq2 approach uses estimated gene counts from 
#       the transcript abundance quantifiers, but not normalized counts.
# Our recommended pipeline for DESeq2 is to use fast transcript abundance 
#       quantifiers upstream of DESeq2
# http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#input-data

# Rscript --vanilla run_deseq2.R $workDir $projectid $organismname \
#         $transcriptome $RscriptsDir "tximport" "kallisto"