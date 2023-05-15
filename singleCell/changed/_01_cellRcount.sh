#!/bin/bash
# -*- ENCODING: UTF-8 -*-

##===============================================================================
## SLURM VARIABLES
#SBATCH --job-name=cellRrun
#SBATCH --cpus-per-task=24
#SBATCH --mem=75Gb
#SBATCH --time=1-00:00:00
#SBATCH -p medium
#SBATCH -o /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.out  
#SBATCH -e /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.err 
##SBATCH --dependency=afterany:571091


# original ram was 258Gb, and cpu 24

# HOW TO RUN ME
#sbatch /home/jmendietaes/programas/pipelines/singleCell/changed/01_cellRcount.sh

# OBJECTIVE
# run cellRanger

#source $CODEBASE/tfcf/setup.sh

#### TO CHANGE
# Main data project file
basepath="/home/jmendietaes/data/2021/singleCell"

# base name of the files to check (separated by space)
# avoid adding _Library.csv, it will be added later
filesCheck="ECCITE8_OP1_d7 ECCITE8_OP1_d9"


# Path where we will store output data
outputPath="${basepath}/allProcessed"

# Path where we have the resource files from cellRanger
resourcePath="${basepath}/additionalFiles/refdata-gex-mm10-2020-A"
rangerPath="${basepath}/allProcessed/rangerFiles"
# path for scripts location
scriptsPath="/home/jmendietaes/programas/PhD"

# id of files were feature files have _ensgs (space separated)
featureExtraId="ECCITE7 ECCITE8"



#### MAIN CODE
# load modules
module load CellRanger/6.1.1

echo ${SLURM_MEM_PER_NODE}
adjustedMem=$(echo "print(int(round($SLURM_MEM_PER_NODE*0.98,0)/1000))" | python3)
# fast check in cases where RAM given in Gb (value will be lower than 1 in most cases)
# 100Gb / 1000 = 0.1
if [ "$adjustedMem" -lt "5" ]; then
    adjustedMem=$(echo "print(int(round($SLURM_MEM_PER_NODE*0.98,0)))" | python3)
fi  

genomePath=$basepath/additionalFiles/omicstmp/newGenome
# If genome doesn't exist, then create it
if [ ! -f "${genomePath}/genome.fa" ]; then
    bash ${scriptsPath}/singleCell/changed/00_PREP_Genome.sh
fi

if [ ! -f "${outputPath}/omicstmp" ]; then
    mkdir -p ${outputPath}/omicstmp
fi


########## BASIC ANALYSIS
cd $outputPath/omicstmp

for id in ${filesCheck}; do

    echo $id
    eccite_id=$(echo $id | cut -d '_' -f 1)
    echo $eccite_id

    # check that feature file exist
    if [ ! -e $rangerPath/${eccite_id}_Features${featureExtra}.csv ]; then
        echo "You forgot to create the features files"
        exit 1
    fi

    # check for files were ensgs should be in features files
    featureExtra=""
    for fe in ${featureExtraId}; do
        if [ "$eccite_id" = "$fe" ] ; then
            featureExtra="_ensgs"
        fi
    done

    echo "Files"
    #cat $rangerPath/${id}_Library.csv

    echo "Features"
    #cat $rangerPath/ECCITE8_Features.csv

    cellranger count --id=$id \
     --no-bam \
     --libraries=$rangerPath/${id}_Library.csv \
     --transcriptome=$basepath/additionalFiles/omicstmp/newGenomeExtended/ \
     --feature-ref=$rangerPath/${eccite_id}_Features${featureExtra}.csv \
     --localcores=${SLURM_CPUS_PER_TASK} \
     --localmem=${adjustedMem} \
     --expect-cells=10000 &> ${id}.log

    mkdir -p ${outputPath}/Data/$id
    mv $id/outs ${outputPath}/Data/$id
done



######### BASIC ANALYSIS without RNA

cd $outputPath/omicstmp

# for id_original in ECCITE7_Lib1Rep1 ECCITE7_Lib1Rep2 ECCITE7_Lib2; do
for id_original in ${filesCheck}; do

    echo $id_original
	id="${id_original}_onlyRNA"
	echo $id

	grep -v "CRISPR Guide Capture" $rangerPath/${id_original}_Library.csv > $rangerPath/${id}_Library.csv
	
	
    echo "Files"
    #cat $rangerPath/${id}_Library.csv
    
	cellranger count --id=$id \
	 --no-bam \
	 --libraries=$rangerPath/${id}_Library.csv \
	 --transcriptome=$basepath/additionalFiles/omicstmp/newGenomeExtended/ \
		 --localcores=${SLURM_CPUS_PER_TASK} \
		 --localmem=${adjustedMem} \
	 --expect-cells=10000 &> ${id}.log
    
	mkdir -p ${outputPath}/Data/onlyRNA/$id
	mv $id/outs ${outputPath}/Data/onlyRNA/$id
done



