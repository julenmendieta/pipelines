#!/bin/bash
# -*- ENCODING: UTF-8 -*-

# HOW TO RUN ME
# bash /home/julen/programas/pipelines/ChIP/04c_chromHMM_stacked.sh

# OBJECTIVE
# Run ChromHMM BinarizeBam on a set of bam files

###############  VARIABLES  ##############
# RAM
memo=30000

# n CPU
nCPU=32

# Tab separated file with chrom name and length
chromosomelengthfile="/home/julen/genomes/mm10_reordered/mm10.reordered.sizes.shortTab"
# Genome assembly id
assembly="mm10"
# ID for the run
runID="stacked"

# Path to directory containing input bam files
inputbamdir="/scratch/julen/ChIP/bamFiles/subsampled/all"

# File indicating cell type, tag, and control
# A tab delimited file where each row contains the cell type, then the 
#  associated mark, then the name of a bam file, and optionally a corresponding 
#  control bam file
cellmarkfiletable="/scratch/julen/ChIP/allData/08_restartProject/ChromHMM/metadata/metada.tsv"

# If you want chromHMM to merge cell and chip id set to "-stacked", otherwise ""
#  stacked models may help differentiate regions with constitutive chromatin 
#  activities from those with cell-type-specific activities.
#  while the stacked model state definitions are more complex, the resulting 
#  genome annotations are simpler and non-overlapping. With the stacked modeling, 
#  each location is simply assigned to one of N universal states, whereas in 
#  the concatenated model, each location is assigned to one of M states in 
#  K cell types
#stacked=""
stacked="-stacked"

# Path to main output dir
outpath="/scratch/julen/ChIP/allData/08_restartProject/ChromHMM/output/stacked"

# The output directory to which the binarized data files should be written. 
# These files will be named CELL_CHROM_binary.txt or 
outputbinarydir="${outpath}/binarize"

# The output directory to which the model data files should be written. 
outputModelDir="${outpath}/model"

# Bin size f the model (default 200bp)
binsize=200

# define number of states to check (array format)
#numstatesAll=( 16 17 18 )
numstatesAll=( 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 )

# This indicates a threshold for the fold enrichment over expected that must 
# be met or exceeded by the observed count in a bin for a present call. This 
#  parameter can be useful when dealing with very deeply and/or unevenly 
#  sequenced data. By default this parameter value is 0
foldthresh=1.5

# Path to ChromHMM jar file
chromJar="/home/julen/programas/ChromHMM_1.23/ChromHMM/ChromHMM.jar"


# Make needed folders
if [ ! -e ${outputbinarydir}/control ]; then
    mkdir -p ${outputbinarydir}/control
    mkdir -p ${outputbinarydir}/signal
fi



#######################  Functions  #####################
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


#######################  RUN  ########################

# Create Step control file if it doesn't exit
stepControl="${outpath}/chromHMMStep_log.txt"
if [ ! -e ${stepControl} ] ; then
    touch ${stepControl}
fi

mkdir -p ${outputModelDir}
####################################
# BinarizeBam – Convert BAM files to binarized data (needed for LearnModel)
####################################
# items in [] are optional
# java -mx${memo}M -jar ChromHMM.jar BinarizeBam [-b binsize][-c controldir]
# [-center][-e offsetend][-f foldthresh][-g signalthresh][-gzip][-n shift]
# [-o outputcontroldir][-p poissonthresh]
# [-paired|-mixed|[-center][-peaks [-i splitrowindex]]]
# [-s offsetstart][-splitcols [-k
# splitcolindex][-m numsplitcols]][-splitrows [-j numsplitbins]][-stacked]
# [-strictthresh][-t outputsignaldir][-u pseudocountcontrol][-w flankwidthcontrol]
# chromosomelengthfile inputbamdir cellmarkfiletable outputbinarydir

echo -e "Starting BinarizeBam -------------------------------------- \n"

# check content of first line of step control file
linec=`sed "1q;d" ${stepControl}`
if [[ ${linec} != "BinarizeBam" ]]; then 
    # To check use in future: foldthresh signalthresh stacked
    java -mx${memo}M -jar ${chromJar} BinarizeBam -b ${binsize} -paired \
        -f ${foldthresh} -gzip ${stacked} \
        -o ${outputbinarydir}/control \
        -t ${outputbinarydir}/signal \
        ${chromosomelengthfile} ${inputbamdir} ${cellmarkfiletable} ${outputbinarydir}

    echo -e "BinarizeBam - done -------------------------------------- \n"
    # store stage control info
    echo "BinarizeBam" > ${stepControl}
else
    echo -e "BinarizeBam - already done before -------------------------------------- \n"
fi


####################################
# LearnModel – Learn chromatin state model from binarized data
####################################
# LearnModel [-b binsize][-color r,g,b][-d convergedelta]
# [-e loadsmoothemission][-f inputfilelist][-gzip][-h informationsmooth]
# [-holdcolumnorder][-holdroworder][-i outfileID][-init information|random|load]
# [-l chromosomelengthfile][-lowmem][-m modelinitialfile][-many][-n numseq][-
# noautoopen][-nobed][-nobrowser][-noenrich][-noimage][-p maxprocessors][-pseudo][-
# printposterior][-printstatebyline][-r maxiterations][-s seed][-scalebeta][-
# splitrows][-stateordering emission|transition][-t loadsmoothtransition] [-u
# coorddir][-v anchorfiledir][-x maxseconds][-z zerotransitionpower]
# inputdir outputdir numstates assembly

echo -e "Starting LearnModel -------------------------------------- \n"

# check content of first line of step control file

for numstates in ${numstatesAll[@]}; do
    # check if the file exists of it was created with a previous bam version 
    fileNotExistOrOlder "${outputModelDir}/${numstates}/webpage_${numstates}_nonStacked.html" "${inputbamdir}/*bam"
    if [[ ${analyse} == "yes" ]]; then
        mkdir -p ${outputModelDir}/${numstates}

        # To check [-holdcolumnorder -printposterior
        # MIGHT BE INTERESTED IN ATHOR MORE ANNOTATION IN ~/programas/ChromHMM_1.23/ChromHMM/COORDS/mm10
        # LIKE REPETITIONS, TAD borders etc.

        java -mx${memo}M -jar ${chromJar} LearnModel -b ${binsize} \
            -gzip -i ${runID} -l ${chromosomelengthfile} -noautoopen -p ${nCPU} \
            ${outputbinarydir} ${outputModelDir}/${numstates} ${numstates} \
            ${assembly} > ${outputModelDir}/${numstates}/leanLog.txt

        echo -e "LearnModel ${numstates} states - done -------------------------------------- \n"
        # store stage control info
    else
        echo -e "LearnModel ${numstates} states - already done before -------------------------------------- \n"
    fi
done