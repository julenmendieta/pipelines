#!/bin/bash
# -*- ENCODING: UTF-8 -*-

##===============================================================================
## SLURM VARIABLES
#SBATCH --job-name=cellRgenome
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=0-02:00:00
#SBATCH -p short
#SBATCH -o /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.out  
#SBATCH -e /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.err 


# HOW TO RUN ME
#sbatch /home/jmendietaes/programas/pipelines/singleCell/changed/00_PREP_Genome.sh

# OBJECTIVE
# merge genome files for cellRanger



# Path where we have the processed files
basepath=/home/jmendietaes/data/2021/singleCell/additionalFiles
# Path where we have the resoruce files from cellRAnger
#resourcePath="/home/jmendietaes/data/2021/singleCell/additionalFiles/refdata-gex-GRCh38-2020-A"
resourcePath="/home/jmendietaes/data/2021/singleCell/additionalFiles/omicstmp/refdata-gex-mm10-2020-A"
genomeId=$(basename ${resourcePath})


# load modules
module load CellRanger/6.1.1

### CREATE GENOME

#source $CODEBASE/tfcf/setup.sh
if [ ! -e $basepath/omicstmp ]; then
	mkdir -p $basepath/omicstmp
fi
cd $basepath/omicstmp

outPath="$basepath/omicstmp/${genomeId}-Extended"
#if [ ! -e $outPath ]; then
#	mkdir -p $outPath
#fi

mkdir -p $resourcePath/extended

# Copy original reference genome and gtf files
cp $resourcePath/fasta/genome.fa $resourcePath/extended/genome.fa
cp $resourcePath/genes/genes.gtf $resourcePath/extended/genes.gtf

# Add GFP and BFP to reference genome and gtf
#x="GFP"
for x in GFP BFP; do
	pathFA=$basepath/refFasta/${x}.fa

	numberBases=$(cat $pathFA | grep -v "^>" | tr -d "\n" | wc -c)
	gtf="$x\tunknown\texon\t1\t${numberBases}\t.\t+\t.\tgene_id xxx${x}xxx; transcript_id xxx${x}xxx; gene_name xxx${x}xxx; gene_biotype xxxprotein_codingxxx;"
	gtf=$(echo $gtf | sed 's/xxx/"/g')

	echo -e $gtf > ${x}.gtf

	cat ${x}.gtf >> $resourcePath/extended/genes.gtf
	fold -w 60 $pathFA >> $resourcePath/extended/genome.fa

	echo -e "" >> $resourcePath/extended/genome.fa

	rm ${x}.gtf
done

grep ">" $resourcePath/extended/genome.fa
tail -5 $resourcePath/extended/genes.gtf
tail -20 $resourcePath/extended/genome.fa

# Create new reference for celRanger
cd $basepath/omicstmp
cellranger mkref --genome="${genomeId}-Extended"  --fasta=$resourcePath/extended/genome.fa --genes=$resourcePath/extended/genes.gtf
