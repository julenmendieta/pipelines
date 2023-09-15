# Script to delete all files from analysis that will be done again
# Add here unique IDs of each sample separated by a space
#bash /home/jmendietaes/programas/pipelines/ChIP/cluster/06_removeFiles_preMerge.sh

removeChip="DM_Hoxa9 DM_Kmt2b DM_Rpb1 Mye_Arid2 Mye_Kmt2b Mye_MPP8 Mye_RING1B \
Mye_Rpb1 Mye_Znhit1"

delete='yes'


basePath=/home/jmendietaes/data/2021/chip/allProcessed
outpath=${basePath}"/furtherAnalysis/08_projectRestart"
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