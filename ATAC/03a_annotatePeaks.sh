#!/bin/bash
# -*- ENCODING: UTF-8 -*-

# HOW TO RUN ME
#bash /home/jmendietaes/programas/pipelines/ChIP/cluster/csaw/02a_annotatePeaks.sh 

# OBJECTIVE
# Annotate peakfiles in a folder

#basePath="/scratch/julen/ChIP/allData/04_subsamplingNoIgG/outdata/csaw"
basePath=$1
#nCPU=16
nCPU=$2

inpath="${basePath}/binnedPeaks"
outpath="${basePath}/Annot/consensus"


# GTF file for annotation (top be consistent with scRNA data)
# Set to FALSE if you wnat HOMER's default UCSC refGene annotation
gtfFile=/scratch/julen/singleCell/cellRanger/genes.gtf
#gtfFile=FALSE

# If we want a column focussed on repeated elements only
# Path to Homer file with repeat element locations
repeatsPath=/home/julen/programas/HOMER/data/genomes/mm10/mm10.repeats
#repeatsPath=FALSE

# path for the location of the pipeline scripts
scriptsPath="/home/julen/programas/PhD"
# species shortcut for MACS
species="mm"
speciesGenome="mm10"

##===============================================================================

# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command failed with exit code $?."' EXIT

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


# get extra parameters for annotation
if [[ $gtfFile == "FALSE" ]]; then
    extraAnnot=""
else
    extraAnnot="-gtf ${gtfFile}"
fi

############################################################
# Annotate consensus peaks with HOMER, and add annotation to boolean output file
############################################################

if [ ! -e ${outpath} ]; then
	mkdir -p ${outpath}
fi

echo -e "Starting consensus peak annotations ------------------------------\n"

consensusFiles=$(find ${inpath}/*tsv -printf "${inpath}/%f\n" | \
            tr '\n' ' ')
chip="allmerged"
for binnedPeaks in ${consensusFiles}; do
    peaktype="binnedPeak"
    prefix=$(basename $binnedPeaks | cut -d '.' -f 1 | sed 's/allChIPCounts_//g')

    prefix="${chip}_${peaktype}_consensusPeaks_${prefix}"

    ## First part
    mergecols=`seq 2 9 | tr '\n' ','`
    expandparam=''



    ## Second part
    # check if the file exists or it was created with a previous peaks version 
    boolAnotMatr=${outpath}/${prefix}.boolean.annotatePeaks.txt
    fileNotExistOrOlder "${boolAnotMatr}" "${binnedPeaks}"
    # this outputs analyse as yes or no in lowercase
    if [[ ${analyse} == "yes" ]]; then

        annotatePeaks.pl \
                ${binnedPeaks} \
                ${speciesGenome} \
                -gid \
                ${extraAnnot} \
                -cpu ${nCPU} \
                -annStats ${outpath}/${prefix}.annotateStats.txt \
                > ${outpath}/${prefix}.annotatePeaks.txt

        cut -f2- ${outpath}/${prefix}.annotatePeaks.txt | \
            awk 'NR==1; NR > 1 {print $0 | "sort -T '.' -k1,1 -k2,2n"}' | \
            cut -f6- > ${outpath}/tmp.txt
        paste ${binnedPeaks} \
            ${outpath}/tmp.txt > ${boolAnotMatr}
        
        # We add a column for annotations regarding repeat elements
        if [[ $repeatsPath != "FALSE" ]]; then
            annotatePeaks.pl \
                    ${binnedPeaks} \
                    ${speciesGenome} \
                    -gid \
                    -ann ${repeatsPath} \
                    -cpu ${nCPU} \
                    > ${outpath}/${prefix}.annotatePeaks_rep.txt

            cut -f2- ${outpath}/${prefix}.annotatePeaks_rep.txt | \
                awk 'NR==1; NR > 1 {print $0 | "sort -T '.' -k1,1 -k2,2n"}' | \
                cut -f7 > ${outpath}/tmp.txt
            rm ${outpath}/${prefix}.annotatePeaks_rep.txt
            # rename header and paste to consensus table
            sed -i 's/Annotation/Repeats Annotation/g' ${outpath}/tmp.txt
            paste ${boolAnotMatr} \
                ${outpath}/tmp.txt > ${outpath}/tmp2.txt
            mv ${outpath}/tmp2.txt ${boolAnotMatr}
        fi
            
    fi
done

echo -e "consensus peak annotations - Finished ---------------------\n"