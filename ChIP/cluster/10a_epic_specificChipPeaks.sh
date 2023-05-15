# Script to call peaks with Epic (ideal for very broad peaks)
# How to run me
# You'll have to ask for extra RAM salloc -c 1 -p short --mem=30000
# bash /home/jmendietaes/programas/pipelines/ChIP/cluster/10a_epic_specificChipPeaks.sh

chips="K9me3"
control="IgG"
species=mm10
bamsPath="/home/jmendietaes/data/2021/chip/allProcessed/bamfiles/valid/08_projectRestart"
outpath=/home/jmendietaes/data/2021/chip/allProcessed/furtherAnalysis/08_projectRestart/peakCalling/Epic2
nCPU=8

cells="DM Mye"
for chip in ${chips}; do
      echo ${chip}
      for cell in ${cells}; do
            echo ${cell}
            control1=$(find ${bamsPath}/${cell}_${control}*bam)
            infile1=$(find ${bamsPath}/${cell}_${chip}*bam)
            epic2 --treatment ${infile1} --control ${control1} --keep-duplicates \
                  --genome ${species} --output ${outpath}/${cell}_${chip}_epic.out \
                  --guess-bampe
            sleep 5
            cat ${outpath}/${cell}_${chip}_epic.out | \
                  cut -f 1,2,3 > ${outpath}/${cell}_${chip}_epic.bed
      done
      echo
done




# outpath=/home/jmendietaes/data/2021/chip/allProcessed/furtherAnalysis/08_projectRestart/peakCalling/Sicer2

# for chip in ${chips}; do
#       echo ${chip}
#       for cell in ${cells}; do
#             echo ${cell}
#             control1=$(find ${bamsPath}/${cell}_${control}*bam)
#             infile1=$(find ${bamsPath}/${cell}_${chip}*bam)
#             sicer --treatment_file ${infile1} --control_file ${control1} \
#                   --species ${species} --output_directory ${outpath} \
#                   --paired_end --cp ${nCPU}
#             sleep 5
#             cat ${outpath}/${cell}_${chip}_epic.out | \
#                   cut -f 1,2,3 > ${outpath}/${cell}_${chip}_epic.bed
#       done
#       echo
# done