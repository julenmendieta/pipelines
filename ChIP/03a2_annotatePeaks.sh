#!/bin/bash
# -*- ENCODING: UTF-8 -*-

# HOW TO RUN ME
#bash /home/jmendietaes/programas/PhD/ChIP/cluster/csaw/03a2_annotatePeaks.sh 

# OBJECTIVE
# Annotate specific peakfile by given gtf file

#basePath="/scratch/julen/ChIP/allData/04_subsamplingNoIgG/outdata/csaw"
binnedPeaks=$1

outpath=$2
#nCPU=16
nCPU=$3

# Genome
speciesGenome=$4
#speciesGenome="mm10"

# GTF file for annotation (top be consistent with scRNA data)
# Set to FALSE if you wnat HOMER's default UCSC refGene annotation
gtfFile=$5
#gtfFile=/scratch/julen/singleCell/cellRanger/mm10-2020-A_genes.gtf
#gtfFile=FALSE

# If we want a column focussed on repeated elements only
# Path to Homer file with repeat element locations
#repeatsPath=/home/julen/programas/HOMER/data/genomes/mm10/mm10.repeats
repeatsPath=$6
# path for the location of the pipeline scripts
scriptsPath="/home/julen/programas/PhD"
# Path were homer genome files are stored
homerGenome="/home/julen/programas/HOMER/.//data/genomes"

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


peaktype="binnedPeak"
prefix=$(basename $binnedPeaks | cut -d '.' -f 1 | sed 's/allChIPCounts_//g')

prefix="${peaktype}_${prefix}"
cell=$( mapLib=(${prefix//_/ }); 
        echo ${mapLib[1]})

# ## First we check if this genome has already been added to homer
# if [[ ! -e ${homerGenome}/${cell} && $speciesGenome == *fasta ]]; then
# 	loadGenome.pl -name ${cell} -fasta ${speciesGenome} -gtf ${gtfFile} -org ${cell}
#     extra2="-organism ${cell}"
# elif [[ $speciesGenome == *fasta ]]; then
#     extra2="-organism ${cell}"
# else
#     extra2=""
# fi


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
            ${extra2} \
            ${extraAnnot} \
            -cpu ${nCPU} \
            -annStats ${outpath}/${prefix}.annotateStats.txt \
            > ${outpath}/${prefix}.annotatePeaks.txt

    cut -f2- ${outpath}/${prefix}.annotatePeaks.txt | \
        awk 'NR==1; NR > 1 {print $0 | "sort -T '.' -k1,1 -k2,2n"}' | \
        cut -f6- > ${outpath}/tmp.txt
    paste ${binnedPeaks} \
        ${outpath}/tmp.txt > ${boolAnotMatr}
    
    # external GTF file annotation gives error in Annotation columns, so I 
    # will get this one appart only if its a default genome
    if [[ $gtfFile != "FALSE" && $speciesGenome != *fasta ]]; then
        annotatePeaks.pl \
                ${binnedPeaks} \
                ${speciesGenome} \
                -gid \
                -cpu ${nCPU} \
                -annStats ${outpath}/${prefix}.annotateStats.txt \
                > ${outpath}/${prefix}.annotatePeaks.txt

        cut -f2- ${outpath}/${prefix}.annotatePeaks.txt | \
            awk 'NR==1; NR > 1 {print $0 | "sort -T '.' -k1,1 -k2,2n"}' | \
            cut -f7 > ${outpath}/tmp.txt
        rm ${outpath}/${prefix}.annotatePeaks.txt
        # rename header and paste to consensus table
        sed -i 's/Annotation/Annotation2/g' ${outpath}/tmp.txt
        paste ${boolAnotMatr} \
            ${outpath}/tmp.txt > ${outpath}/tmp2.txt
        mv ${outpath}/tmp2.txt ${boolAnotMatr}
    fi

    # We add a column for annotations regarding repeat elements
    if [[ $repeatsPath != "FALSE" ]]; then
        annotatePeaks.pl \
                ${binnedPeaks} \
                ${speciesGenome} \
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


echo -e "consensus peak annotations - Finished ---------------------\n"