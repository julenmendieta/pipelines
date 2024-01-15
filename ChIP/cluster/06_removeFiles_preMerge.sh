# Script to delete all files from analysis that will be done again
# Add here unique IDs of each sample separated by a space
#bash /home/jmendietaes/programas/pipelines/ChIP/cluster/06_removeFiles_preMerge.sh

removeChip="Mye-HD5_Brd9_ChIP19 Mye-HD5_H3K4me3_ChIP19 Mye-HD5_Kmt2d_ChIP19 \
Mye_Brd9_ChIP19 Mye_H3K4me3_ChIP19 Mye_Kmt2d_ChIP19"

delete='yes'


basePath=/home/jmendietaes/data/2021/chip/allProcessed
outpath=${basePath}"/furtherAnalysis/14_hexadienol"
if [ ! -e ${outpath}/peakCalling/MACS2/peaks/mergedReplicates ]; then
    mkdir -p ${outpath}/peakCalling/MACS2/peaks/mergedReplicates
    mkdir -p ${outpath}/HOMER/peakAnnotation/mergedReplicates
    #mkdir -p ${outpath}/peakCalling/MACS2/consensusPeaks/bySameChip/mergedReplicates
fi

for re in ${removeChip}; do
    ls ${outpath}/peakCalling/MACS2/peaks/${re}* | tr ' ' '\n'
    ls ${outpath}/HOMER/peakAnnotation/${re}* | tr ' ' '\n'
    if [[ ${delete} == "yes" ]]; then
        rm ${outpath}/peakCalling/MACS2/peaks/${re}*
        #rm ${outpath}/peakCalling/MACS2/consensusPeaks/bySameChip/${re}*
        rm ${outpath}/HOMER/peakAnnotation/${re}*
    else
        mv ${outpath}/peakCalling/MACS2/peaks/${re}* ${outpath}/peakCalling/MACS2/peaks/mergedReplicates
        #mv ${outpath}/peakCalling/MACS2/consensusPeaks/bySameChip/${re}* ${outpath}/peakCalling/MACS2/consensusPeaks/bySameChip/mergedReplicates
        mv ${outpath}/HOMER/peakAnnotation/${re}* ${outpath}/HOMER/peakAnnotation/mergedReplicates
    fi
done