# this was a fast test to get CPM in selected regions of interest and then compare them


bamfiles="/home/jmendietaes/data/2021/chip/allProcessed/bamfiles/valid/subsampled/DM_Smarcb1-merged-sub173331792.sort.rmdup.rmblackls.rmchr.bam \
/home/jmendietaes/data/2021/chip/allProcessed/bamfiles/valid/subsampled/DM_Kmt2a-merged-sub101786848.sort.rmdup.rmblackls.rmchr.bam \
/home/jmendietaes/data/2021/chip/allProcessed/bamfiles/valid/subsampled/DM_Kmt2d-merged-sub146746296.sort.rmdup.rmblackls.rmchr.bam \
/home/jmendietaes/data/2021/chip/allProcessed/bamfiles/valid/DM_Stat5a_Cyn1_S8.sort.rmdup.rmblackls.rmchr.bam \
/home/jmendietaes/data/2021/chip/allProcessed/bamfiles/valid/subsampled/Mye_Smarcb1-merged-sub173331792.sort.rmdup.rmblackls.rmchr.bam \
/home/jmendietaes/data/2021/chip/allProcessed/bamfiles/valid/subsampled/Mye_Kmt2a-merged-sub101786848.sort.rmdup.rmblackls.rmchr.bam \
/home/jmendietaes/data/2021/chip/allProcessed/bamfiles/valid/subsampled/Mye_Kmt2d-merged-sub146746296.sort.rmdup.rmblackls.rmchr.bam \
/home/jmendietaes/data/2021/chip/allProcessed/bamfiles/valid/Mye_Stat5a_Cyn1_S7.sort.rmdup.rmblackls.rmchr.bam"
outpath=/home/jmendietaes/data/2021/chip/allProcessed/furtherAnalysis/07_leukPaper/addStat5
featureCounts \
        -F SAF \
        -O \
        --fracOverlap 0.2 \
        -T 4 \
        -p --donotsort \
        -a ${outpath}/test.saf \
        -o ${outpath}/stat5/stat5_narrowPeak_consensusPeaks.featureCounts.txt \
        ${bamfiles}

prevFields=6
lengthCol=5
scriptsPath="/home/jmendietaes/programas/PhD"
export PATH="/home/jmendietaes/programas/miniconda3/envs/Renv/bin:$PATH"
export PATH="/home/jmendietaes/programas/miniconda3/bin:$PATH"

# set the output file path and copy the content of original
featureCPM=${outpath}/stat5/stat5_narrowPeak_consensusPeaks.featureCounts.CPM.txt

echo -e "FileName\tSampleName\tCellType\tStatus" > ${outpath}/sampleInfo.txt
for bam in ${bamfiles}; do
    sname=`basename $bam  | sed 's/.sort.rmdup.rmblackls.rmchr.bam//g'`
    #sname=$bam
    cell=(${sname//_/ }) ; cell=${cell[0]}
    status=(${sname//_/ }) ; status=${status[1]}; status=(${status//-/ }); status=${status[0]}
    echo -e ${bam}"\t"${sname}"\t"${cell}"\t"${status} >> ${outpath}/sampleInfo.txt
done

Rscript ${scriptsPath}/ChIP/cluster/02_NR_metric_CPMfromFeatureCounts.r \
                        -i ${outpath}/stat5/stat5_narrowPeak_consensusPeaks.featureCounts.txt \
                        -o ${featureCPM} \
                        -s ${outpath}/sampleInfo.txt \
                        -l ${lengthCol} \
                        -d ${prevFields}

python ${scriptsPath}/ChIP/cluster/02_NR_gatherAllInTable.py \
            ${outpath}/stat5 \
            ${outpath}/stat5 \
            ${outpath}/stat5


bamfiles="/home/jmendietaes/data/2021/chip/allProcessed/bamfiles/valid/GMPvitro_Kmt2a_200721_S10.sort.rmdup.rmblackls.rmchr.bam /home/jmendietaes/data/2021/chip/allProcessed/bamfiles/valid/LSKvitro_Kmt2a_170821_S31.sort.rmdup.rmblackls.rmchr.bam /home/jmendietaes/data/2021/chip/allProcessed/bamfiles/valid/mergedReplicates/MEPvitro_Kmt2a_170821_S30.sort.rmdup.rmblackls.rmchr.bam /home/jmendietaes/data/2021/chip/allProcessed/bamfiles/valid/mergedReplicates/MEPvitro_Kmt2a_200721_S11.sort.rmdup.rmblackls.rmchr.bam"
outpath=/home/jmendietaes/data/2021/chip/allProcessed/furtherAnalysis/lowInputCorrelation/counting
featureCounts \
        -F SAF \
        -O \
        --fracOverlap 0.2 \
        -T 4 \
        -p --donotsort \
        -a ${outpath}/Kmt2a.saf \
        -o ${outpath}/Kmt2a/Kmt2a.featureCounts.txt \
        ${bamfiles}

prevFields=6
lengthCol=5
scriptsPath="/home/jmendietaes/programas/PhD"
export PATH="/home/jmendietaes/programas/miniconda3/envs/Renv/bin:$PATH"

# set the output file path and copy the content of original
featureCPM=${outpath}/Kmt2a/Kmt2a.featureCounts.CPM.txt

echo -e "FileName\tSampleName\tCellType\tStatus" > ${outpath}/sampleInfo.txt
for bam in ${bamfiles}; do
    sname=`basename $bam  | sed 's/.sort.rmdup.rmblackls.rmchr.bam//g'`
    #sname=$bam
    cell=(${sname//_/ }) ; cell=${cell[0]}
    status=(${sname//_/ }) ; status=${status[1]}; status=(${status//-/ }); status=${status[0]}
    echo -e ${bam}"\t"${sname}"\t"${cell}"\t"${status} >> ${outpath}/sampleInfo.txt
done

Rscript ${scriptsPath}/ChIP/cluster/02_NR_metric_CPMfromFeatureCounts.r \
                        -i ${outpath}/Kmt2a/Kmt2a.featureCounts.txt \
                        -o ${featureCPM} \
                        -s ${outpath}/sampleInfo.txt \
                        -l ${lengthCol} \
                        -d ${prevFields}


                    


bamfiles="/home/jmendietaes/data/2021/chip/allProcessed/bamfiles/valid/mergedReplicates/GMPvitro_Kmt2d_170821_S50.sort.rmdup.rmblackls.rmchr.bam /home/jmendietaes/data/2021/chip/allProcessed/bamfiles/valid/mergedReplicates/GMPvitro_Kmt2d_200721_S5.sort.rmdup.rmblackls.rmchr.bam /home/jmendietaes/data/2021/chip/allProcessed/bamfiles/valid/mergedReplicates/LSKvitro_Kmt2d_170821_S49.sort.rmdup.rmblackls.rmchr.bam /home/jmendietaes/data/2021/chip/allProcessed/bamfiles/valid/mergedReplicates/LSKvitro_Kmt2d_200721_S7.sort.rmdup.rmblackls.rmchr.bam /home/jmendietaes/data/2021/chip/allProcessed/bamfiles/valid/MEPvitro_Kmt2d_200721_S6.sort.rmdup.rmblackls.rmchr.bam"
outpath=/home/jmendietaes/data/2021/chip/allProcessed/furtherAnalysis/lowInputCorrelation/counting
featureCounts \
        -F SAF \
        -O \
        --fracOverlap 0.2 \
        -T 4 \
        -p --donotsort \
        -a ${outpath}/Kmt2d.saf \
        -o ${outpath}/Kmt2d/Kmt2d.featureCounts.txt \
        ${bamfiles}

prevFields=6
lengthCol=5
scriptsPath="/home/jmendietaes/programas/PhD"
export PATH="/home/jmendietaes/programas/miniconda3/envs/Renv/bin:$PATH"

# set the output file path and copy the content of original
featureCPM=${outpath}/Kmt2d/Kmt2d.featureCounts.CPM.txt

echo -e "FileName\tSampleName\tCellType\tStatus" > ${outpath}/sampleInfo.txt
for bam in ${bamfiles}; do
    sname=`basename $bam  | sed 's/.sort.rmdup.rmblackls.rmchr.bam//g'`
    #sname=$bam
    cell=(${sname//_/ }) ; cell=${cell[0]}
    status=(${sname//_/ }) ; status=${status[1]}; status=(${status//-/ }); status=${status[0]}
    echo -e ${bam}"\t"${sname}"\t"${cell}"\t"${status} >> ${outpath}/sampleInfo.txt
done

Rscript ${scriptsPath}/ChIP/cluster/02_NR_metric_CPMfromFeatureCounts.r \
                        -i ${outpath}/Kmt2d/Kmt2d.featureCounts.txt \
                        -o ${featureCPM} \
                        -s ${outpath}/sampleInfo.txt \
                        -l ${lengthCol} \
                        -d ${prevFields}


                    