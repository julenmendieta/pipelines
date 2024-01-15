#!/bin/bash
# -*- ENCODING: UTF-8 -*-

##==============================================================================
## SLURM VARIABLES
#SBATCH --job-name=TEcopies
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --time=00-24:00:00
#SBATCH -p short
#SBATCH -o /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.out  
#SBATCH -e /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.err 
##SBATCH --dependency=afterany:793158


##SBATCH --mail-type=END
##SBATCH --mail-user=user@mail.es
# HOW TO RUN ME
# sbatch /home/jmendietaes/programas/pipelines/mRNA/cluster/03c_TEcopies.sh \
#/home/jmendietaes/data/2021/mRNA/allProcessed/TEtranscripts \
#03_TEtranscripts

# OBJECTIVE
# Count transcripts on TE copies

# NOTES
# Code based on notes from https://doi.org/10.1038/s41556-021-00785-9


##==============================================================================
## GLOBAL VARIABLES


## path to the main fully processed data output
basePath=$1
#basePath="/home/jmendietaes/data/2021/mRNA/allProcessed/TEtranscripts"

## Extra path varaibles
# GTF to TE
TEgtf=/home/jmendietaes/referenceGenomes/annotRegions/mm10/mm10_rmsk_TE.gtf

## Project name
# We will find salmon output and STAR BAM files in there using our defined
#   folder structure
projectName=$2
#projectName=03_TEtranscripts

# Coma separated string with all posible control IDs 
# Used for batch corrected analysis of same KOs
# Set to "no" to not do it
posibleControls="NTC"
#posibleControls="no"


# Path to extra scripts
extraScriptsPath="/home/jmendietaes/programas/pipelines/mRNA/cluster/sub-scripts"

# Where to look for Salmon count files
# A path with folder named by sample and containing quant.sf inside
# I link the interest folders from allProcessed/counts/salmon to
#       another location to divide data by projects
#salmonC="${basePath}/countGroups/${projectName}"
# Were to look for bam files
bamsPath="${basePath}/bamfiles/wholeGenome/valid/${projectName}"
outpath=${basePath}"/furtherAnalysis/${projectName}"



##==============================================================================
## Required Software
module load SAMtools/1.12-GCC-10.2.0
# GCC needed for STAR
module load GCCcore/11.2.0
# this is for bamCoverage, but i already have it in conda and newer version
#module load deepTools/3.2.0-foss-2018b-Python-2.7.15


# Path to STAR binary
TEtranscriptsP="/home/jmendietaes/programas/miniconda3/bin/TEtranscripts"
TEcountP="/home/jmendietaes/programas/miniconda3/bin/TEcount"
# Export R path
export PATH="/home/jmendietaes/programas/miniconda3/envs/Renv/bin:$PATH"

module load HTSeq/0.11.0-foss-2018b-Python-2.7.15
module load Python/2.7.15-foss-2018b


# Variables from the job
nCPU=$SLURM_CPUS_PER_TASK
##==============================================================================

# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

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

##==============================================================================
# We filtered for uniquely mapping RNA-seq reads and unravelled specific 
#  TE copies with expression changes after 5-FU treatment at H16, D3 and D10. 
#  This analysis should be taken with caution, as many recent and potentially 
#  active TE copies will not be included due to mapping issues
# Multimapped reads were filtered with Samtools MAPQ > 50 to extract the 
#  uniquely mapped reads. The HT-seq count (v.0.5.4p3.) algorithm62 was applied 
#  to assign aligned reads to the genomic instances of TE copies using the 
#  following command line ‘htseq-count --s no --m intersection --nonempty’. 
#  Annotation files were constructed from RepeatMasker 
#  (http://www.repeatmasker.org). Differentially expressed TE copies were 
#  identified with the use of the DESeq R package


if [ ! -e ${outpath}/TEtranscript ]; then
    mkdir -p ${outpath}/TEtranscript
fi

#################################
## Define strandness of samples
#################################

allbams=$(find ${bamsPath}/*bam -printf "${bamsPath}/%f\n")

samples=$(for i in $bamsPath/*bam; do echo $i | sed "s+${bamsPath}\/++g" | \
                                       sed 's/.STAR.sort.rmchr.bam//g'; done)

lTypes=$(for s in ${samples}; do 
    grep -A 1 "LIBRARY TYPE" ${basePath}/QC/allStepsSummary/summary_${s}_TE.txt |\
    tail -n 1; done | sort | uniq)
nTypes=$(echo ${lTypes} | wc -w)
if [ "$nTypes" -gt 1 ] ; then
    echo "ERROR: More than one library type"
    echo $nTypes 
    exit 1
fi
# Define strandness
# 0 (unstranded), 1 (stranded) and 2 (reversely stranded)
lType_="no"
if [[ ${lTypes} == *"SF" ]] ; then
    lType_="yes"
elif [[ ${lTypes} == *"SR" ]] ; then
    lType_="reverse"
fi


#################################
## Count TEs
#################################

#samtools view -hb -q 50 --threads ${nCPU} ${bam} > filtered.bam

# Check if we have already a log for this check
fileNotExistOrOlder "${outpath}/TEtranscript/${projectName}_TEcopies.txt" \
                "${allbams} ${TEgtf}"
if [[ ${analyse} == "yes" ]]; then
    # Add header info
    header=$(echo "TEcopy ${samples}" | tr '\n' ' ' | \
            tr ' ' '\t')
    echo -e "${header}" > ${outpath}/TEtranscript/${projectName}_TEcopies.txt
    # We filter by MAPQ >= 50 directly in htseq-count
    htseq-count ${allbams} ${TEgtf} --format "bam" --order "pos" --stranded ${lType_} \
        --minaqual 50 --idattr transcript_id -m intersection-nonempty >> \
        ${outpath}/TEtranscript/${projectName}_TEcopies.txt
else
    echo "htseq-count already ran before"

fi

#################################
## DESeq2
#################################

# Check if we have already a log for this check
fileNotExistOrOlder "${outpath}/TEtranscript/${projectName}_TEcopies_DESeq2.txt" \
                "${outpath}/TEtranscript/${projectName}_TEcopies.txt"
if [[ ${analyse} == "yes" ]]; then

    # Get FC values
    Rscript ${extraScriptsPath}/03a_NR_TEcopiesDeseq2.r \
            --input_file ${outpath}/TEtranscript/${projectName}_TEcopies.txt \
            --outfile ${outpath}/TEtranscript/${projectName}_TEcopies_DESeq2.txt \
            --cores ${nCPU} \
            --controls ${posibleControls}
else
    echo "DESeq2 already ran before"

fi


echo -e "END ------------------------------------------------------------------"

seff $SLURM_JOBID

exit 0