# Script to delete all files from analysis that will be done again
# Add here unique IDs of each sample separated by a space
# bash /home/jmendietaes/programas/pipelines/ATAC-KO/cluster/06_removeFiles_preMerge.sh
removeChip="DMd4-*_ATAC22_"

delete='yes'

basePath=/home/jmendietaes/data/2021/ATAC/allProcessed
outpath=${basePath}"/furtherAnalysis/08e_koATAC"
if [ ! -e ${outpath}/peakCalling/MACS2/peaks/mergedReplicates ]; then
    mkdir -p ${outpath}/peakCalling/MACS2/peaks/mergedReplicates
    mkdir -p ${outpath}/HOMER/peakAnnotation/mergedReplicates
    mkdir -p ${outpath}/peakCalling/MACS2/consensusPeaks/bySameChip/mergedReplicates
fi

for re in ${removeChip}; do
    ls ${outpath}/peakCalling/MACS2/peaks/${re}* | tr ' ' '\n'
    ls ${outpath}/HOMER/peakAnnotation/${re}* | tr ' ' '\n'
    if [[ ${delete} == "yes" ]]; then
        rm ${outpath}/peakCalling/MACS2/peaks/${re}*
        rm ${outpath}/peakCalling/MACS2/consensusPeaks/bySameChip/${re}*
        rm ${outpath}/HOMER/peakAnnotation/${re}*
    else
        mv ${outpath}/peakCalling/MACS2/peaks/${re}* ${outpath}/peakCalling/MACS2/peaks/mergedReplicates
        mv ${outpath}/peakCalling/MACS2/consensusPeaks/bySameChip/${re}* ${outpath}/peakCalling/MACS2/consensusPeaks/bySameChip/mergedReplicates
        mv ${outpath}/HOMER/peakAnnotation/${re}* ${outpath}/HOMER/peakAnnotation/mergedReplicates
    fi
done