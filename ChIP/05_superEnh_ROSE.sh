#!/bin/bash
# -*- ENCODING: UTF-8 -*-

# Following steps at 
# https://www.sciencedirect.com/science/article/pii/S221501612030385X

# HOW TO RUN ME
# bash ...sh

# This scripts runs ROSE to get Super Enhancer annotation from K27ac data

PATHTO=/home/julen/programas/ROSE
PYTHONPATH=$PATHTO/lib
export PYTHONPATH
export PATH=$PATH:$PATHTO/bin

ROSEdir="/home/julen/programas/ROSE"
k27Peaks="/scratch/julen/ChIP/allData/08_restartProject/peaks/Mye_K27ac-merged_peaks.broadPeak"
k27Bam="/scratch/julen/ChIP/bamFiles/notSubsampled/Mye/Mye_K27ac-merged.sort.rmdup.rmblackls.rmchr.bam"
controlBam="/scratch/julen/ChIP/bamFiles/IgG/Mye_IgG-merged.sort.rmdup.rmblackls.rmchr.bam"
outpath="/scratch/julen/ChIP/allData/08_restartProject/outdata/other/superEnh/Mye"
ID="Mye"
# Maximum linking distance for stitching
STITCH=12500
# Transcription Start Size Window
TSS=2000
ROSE=${ROSEdir}/bin/ROSE_main.py


cd ${ROSEdir}
peakGFF=${outpath}/${ID}_ROSE.txt
awk 'BEGIN{OFS="\t"} {print $1, $4, "peak", $2, $3, $5, $6, ".", "ID="NR}' ${k27Peaks} > ${peakGFF}
python ${ROSE} -g mm10 -i ${peakGFF} -r ${k27Bam} -c ${controlBam} \
    -o ${outpath} -s ${STITCH} -t ${TSS}



# arser.add_option("--custom", dest="custom_genome", default=None,
#         help = "Enter the custom genome annotation refseq.ucsc")

# #optional flags
# parser.add_option("-b","--bams", dest="bams",nargs = 1, default=None,
#         help = "Enter a comma separated list of additional bam files to map to")


# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3841062/
# ROSE was run with a
# stitching distance of 12,500 bp, i.e. we allowed enhancers within 12,500 bp to
# be stitched together. In addition, we used a promoter exclusion zone of 2,000
# bp, i.e. if a constituent enhancer was wholly contained within a window +/-
# 1,000 bp around an annotated transcription start site, the constituent enhancer
# was excluded from stitching.
# Stitched enhancers were assigned to the expressed transcript whose TSS was
# the nearest to the center of the stitched enhancer. Expressed transcripts were
# defined as having an at least 0.5 mean rpm/bp H3K27ac ChIP-Seq density in a
# window 500 bases up- and downstream of the TSS.

# Thus, the only defining feature of super-enhancers is an exceptionally 
# high degree of enrichment of transcriptional activators or chromatin marks 
# as determined by ChIP-seq, which is assessed in step 3.

#  Super-enhancers in multiple myeloma cells and other tumors were strongly 
#  enriched for binding of BRD4 (refs. 24,29), a coactivator that can bind 
#  acetylated histones and that interacts directly with the Mediator complex 
#  and elongation factors30,31. 