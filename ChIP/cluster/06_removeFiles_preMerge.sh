# Script to delete all files from analysis that will be done again
# Add here unique IDs of each sample separated by a space
removeChip="MEPvitro_Men1_170821"
basePath=/home/jmendietaes/data/2021/chip/allProcessed
outpath=${basePath}"/furtherAnalysis/subsampled_noIgG"
if [ ! -e ${outpath}/peakCalli
ng/MACS2/peaks/mergedReplicates ]; then
    mkdir -p ${outpath}/peakCalling/MACS2/peaks/mergedReplicates
    mkdir -p ${outpath}/HOMER/peakAnnotation/mergedReplicates
fi

for re in ${removeChip}; do
    ls ${outpath}/peakCalling/MACS2/peaks/${re}* | tr ' ' '\n'
    mv ${outpath}/peakCalling/MACS2/peaks/${re}* ${outpath}/peakCalling/MACS2/peaks/mergedReplicates
    ls ${outpath}/HOMER/peakAnnotation/${re}* | tr ' ' '\n'
    mv ${outpath}/HOMER/peakAnnotation/${re}* ${outpath}/HOMER/peakAnnotation/mergedReplicates
done