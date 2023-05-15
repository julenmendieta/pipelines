#!/bin/bash
# -*- ENCODING: UTF-8 -*-

##===============================================================================
## SLURM VARIABLES
#SBATCH --job-name=hichip_loop
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=05:00:00
#SBATCH -p short
#SBATCH -o /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.out  
#SBATCH -e /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.err 


# HOW TO RUN ME
# sbatch /home/jmendietaes/programas/pipelines/HiChIP/cluster/03a_Dovetail_loopCalling.sh \
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
baseConfig=/home/jmendietaes/data/2021/HiChIP/configFiles
binSize=2500
FiHiChIP="/home/jmendietaes/programas/FitHiChIP/FitHiChIP_HiCPro.sh"

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
outpair="${basePath}/pairtools"

## load modules
hichipScripts="/home/jmendietaes/programas/pipelines/HiChIP/cluster"
export PATH="/home/jmendietaes/programas/miniconda3/envs/Renv/bin:$PATH"
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

if [ ! -e ${outpair}/preLoop ]; then
    mkdir -p ${outpair}/preLoop
fi
if [ ! -e "${basePath}/outLoop" ]; then
    mkdir -p "${basePath}/outLoop"
fi




# Get list of all valid pairfiles
allPairs=$(find ${outpair}/*gz -printf "${outpair}/%f\n" | \
            tr '\n' ' ')
allLabels=`for i in $allPairs; do basename ${i} | cut -d '.' -f 1 | cut -d '_' -f 1,2,3; done | tr '\n' ' '`

for peaktype in broadPeak; do 
    echo ${peaktype}
    for ap in ${allPairs}; do
        echo ${ap}
        label=$(echo $(basename ${ap}) | cut -d '.' -f 1 | cut -d '_' -f 1,2)
        peakfile=$(find ${chipPath}/${label}*${peaktype} )

        #########################
        ## Format pairs file
        #########################
        echo "Format pairs file"
        # check if the file exists of it was created with a previous bam version 
        fileNotExistOrOlder "${outpair}/preLoop/${label}.pairs.gz" \
                            "${ap}"
       
        # this outputs analyse as yes or no in lowercase
        if [[ ${analyse} == "yes" ]]; then
            zcat ${ap} | grep -v '#' | \
                awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5"\t"$7}' | \
                gzip -c > ${outpair}/preLoop/${label}.pairs.gz
        fi

        #########################
        ## Modify config file and run FiHiChIP
        #########################
        echo "Modify config file and run FiHiChIP"
        # check if the file exists of it was created with a previous bam version 
        
        ## First of all we have to modify the example config file
        configf=${baseConfig}/modified/${label}_config.txt 

        fileNotExistOrOlder "${baseConfig}/modified/${label}_config.txt" \
                        "${outpair}/preLoop/${label}.pairs.gz"

        if [[ ${analyse} == "yes" ]]; then
            cp ${baseConfig}/example_config.txt ${configf}

            # Change file indicating input .pairs file path
            sed -i "s|ValidPairs_cc|${outpair}/preLoop/${label}\.pairs\.gz|g" ${configf}
            # Change line containing input ChIP coordinate file path
            sed -i "s|PeakFile_cc|${peakfile}|g" ${configf}
            # Change line containing output directory name
            sed -i "s|OutDir_cc|${basePath}/outLoop/${label}_${binSize}bp|g" ${configf}
            # Change bin size 
            sed -i "s/BINSIZE_cc/${binSize}/g" ${configf}
            # Change reference genome chromosome sizes
            sed -i "s|ChrSizeFile_cc|${REFERENCE_DIR}\.sizes\.shortTab|g" ${configf}
            # Change output prefix
            sed -i "s/PREFIX_cc/${label}_${binSize}bp/g" ${configf}
        fi

        fileNotExistOrOlder "${basePath}/outLoop/${label}_${binSize}bp/Summary_results_FitHiChIP.html" \
                            "${outpair}/preLoop/${label}.pairs.gz"

        if [[ ${analyse} == "yes" ]]; then
            # Run FiHiChIP
            ${FiHiChIP} -C ${configf}
        fi

    done
done

