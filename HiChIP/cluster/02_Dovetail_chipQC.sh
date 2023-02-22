#!/bin/bash
# -*- ENCODING: UTF-8 -*-

##===============================================================================
## SLURM VARIABLES
#SBATCH --job-name=hichip_qc
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=05:00:00
#SBATCH -p short
#SBATCH -o /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.out  
#SBATCH -e /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.err 


# HOW TO RUN ME
# sbatch /home/jmendietaes/programas/PhD/HiChIP/cluster/02_Dovetail_chipQC.sh \
#/home/jmendietaes/data/2021/HiChIP/allProcessed \
#/home/jmendietaes/referenceGenomes/mm10_reordered/mm10.reordered
#/home/jmendietaes/data/2021/chip/allProcessed/furtherAnalysis/08_projectRestart/peakCalling/MACS2/peaks

# OBJECTIVE
# QC of HiChIP using ChIP data as stated in https://hichip.readthedocs.io/en/latest

## GLOBAL VARIABLES
## path to the main fully processed data output
basePath=$1
#basePath="/home/jmendietaes/data/2021/HiChIP/allProcessed"
# here we will stored the final filtered bam files
bamsPath="${basePath}/bamfiles/valid"

##  REFERENCE Genome info
REFERENCE_DIR=$2
#REFERENCE_DIR="/home/jmendietaes/referenceGenomes/mm10_reordered/mm10.reordered"
# REFERENCE_DIR includes the path and naming base of the reference genome, we will
#   add .sizes, .blacklist.bed to the base to get the rest of paths. Genome
#   indexes must follow the same base name


## Path to MACS output of ChIP-seq experiments
chipPath=$3
#chipPath="/home/jmendietaes/data/2021/chip/allProcessed/furtherAnalysis/\
#08_projectRestart/peakCalling/MACS2/peaks"

## Other variables based on the previous
GenomeIndex=$REFERENCE_DIR
statsOut="${basePath}/chipIntersection"
#chrOrder=$REFERENCE_DIR"/mm10_Bowtie2/names.txt"
#chr_genome_size=$REFERENCE_DIR".sizes"
#BlackList=$REFERENCE_DIR".blacklist.bed"
#wigToBigWig="/home/jmendietaes/programas/binPath/wigToBigWig"
#picardPath='/home/jmendietaes/programas/picard/picard.jar'

## load modules
hichipScripts="/home/jmendietaes/programas/PhD/HiChIP/cluster"
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

if [ ! -e ${statsOut} ]; then
    mkdir -p ${statsOut}
fi


# Get list of all valid bamfiles
allbams=$(find ${bamsPath}/*bam -printf "${bamsPath}/%f\n" | \
            tr '\n' ' ')
allLabels=`for i in $allbams; do basename ${i} | cut -d '.' -f 1 | cut -d '_' -f 1,2,3; done | tr '\n' ' '`
# Get short labels to look for ChIP experiments
allLabels_short=$(\
    for filename in ${allLabels}; do 
        mapLib=(${filename//_/ }); 
        echo ${mapLib[0]}_${mapLib[1]}; done )

for peaktype in narrowPeak broadPeak; do 
    for label in ${allLabels_short}; do
        peakfile=$(find ${chipPath}/${label}*${peaktype} )
        bamfile=$(echo $allbams | tr ' ' '\n' | grep ${label})

        # check if the file exists of it was created with a previous bam version 
        fileNotExistOrOlder "${statsOut}/${label}_${peaktype}_hichip_qc_metrics.txt" \
                            "${peakfile} ${bamfile}"
        # Get overlap statistics
        #Metric	                                Shallow Seq (20M)	Deep Seq (100-200M)
        #Tot reads 1000 bp around peak cntr	    >2%	            >2%
        # this outputs analyse as yes or no in lowercase
        if [[ ${analyse} == "yes" ]]; then
            ${hichipScripts}/01_NR_enrichment_stats.sh -g ${REFERENCE_DIR}.fa \
                            -b ${bamfile} -p ${peakfile} -t ${SLURM_CPUS_PER_TASK} \
                            -x ${statsOut}/${label}_${peaktype}
        fi


        # Plot overlap
        fileNotExistOrOlder "${statsOut}/${label}_${peaktype}_enrichment.png" \
                            "${peakfile} ${bamfile}"
        if [[ ${analyse} == "yes" ]]; then

            python3 ${hichipScripts}/01_NR_plot_chip_enrichment.py -bam ${bamfile} \
                        -peaks ${peakfile} \
                        -output ${statsOut}/${label}_${peaktype}_enrichment.png
            

        fi
    done
done

# If requirements passed
#Metric	                                Shallow Seq (20M)	Deep Seq (100-200M)
#No-Dup Read Pairs	                    >75%	            >50%
#No-dup cis read pairs ≥ 1kb	        >20%	            >20%
#Tot reads 1000 bp around peak cntr	    >2%	                >2%

#For shallow sequenced libraries
#   - proceed to deep sequencing (~150 M read pairs per library) 
#For deep sequencing 
#   – proceed with downstream analyses



# Convert filtered pairs file to Hi-C Pro valid pairs format
grep -v '#' <*.pairs> | \
    awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5"\t"$7}' | \
    gzip -c > <output.pairs.gz>

