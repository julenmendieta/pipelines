# Script to delete all files from analysis that will be done again
# Add here unique IDs of each sample separated by a space
removeChip="Activated-Brd9_ATAC2 Activated-Egr2_ATAC2 \
Activated-NtC5_ATAC2 Quiescent-Brd9_ATAC2 Quiescent-Egr2_ATAC2 \
Quiescent-NtC5_ATAC2"
basePath=/home/jmendietaes/data/2021/ATAC/allProcessed
outpath=${basePath}"/furtherAnalysis/05_laura"
if [ ! -e ${outpath}/peakCalling/MACS2/peaks/mergedReplicates ]; then
    mkdir -p ${outpath}/peakCalling/MACS2/peaks/mergedReplicates
    mkdir -p ${outpath}/HOMER/peakAnnotation/mergedReplicates
fi

for re in ${removeChip}; do
    ls ${outpath}/peakCalling/MACS2/peaks/${re}* | tr ' ' '\n'
    mv ${outpath}/peakCalling/MACS2/peaks/${re}* ${outpath}/peakCalling/MACS2/peaks/mergedReplicates
    ls ${outpath}/HOMER/peakAnnotation/${re}* | tr ' ' '\n'
    mv ${outpath}/HOMER/peakAnnotation/${re}* ${outpath}/HOMER/peakAnnotation/mergedReplicates
done