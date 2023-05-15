#!/bin/bash
# -*- ENCODING: UTF-8 -*-

##===============================================================================
## SLURM VARIABLES
#SBATCH --job-name=peaksVScoverage
#SBATCH --cpus-per-task=16
#SBATCH --mem=15G
#SBATCH --time=02-10:00:00
#SBATCH -p medium
#SBATCH -o /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.out  
#SBATCH -e /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.err 


# HOW TO RUN ME
# sbatch /home/jmendietaes/programas/pipelines/ChIP/cluster/03_peaksVScoverage.sh 

# load modules
module load Sambamba/0.7.0
module load SAMtools/1.12-GCC-10.2.0
module load MACS2/2.2.7.1-foss-2018b-Python-3.6.6
module load BEDTools/2.27.1-foss-2018b

##### TO CHANGE #####
species="mm"
inPath=""
outPath="/home/jmendietaes/data/2021/chip/allProcessed/furtherAnalysis"

myeControl="/home/jmendietaes/data/2021/chip/allProcessed/bamfiles/valid/Mye_IgG_31032021_S22.sort.rmdup.rmblackls.rmchr.bam"
dmControl="/home/jmendietaes/data/2021/chip/allProcessed/bamfiles/valid/DM_input_140421_S36.sort.rmdup.rmblackls.rmchr.bam"
chipControl="/home/jmendietaes/data/2021/chip/allProcessed/bamfiles/valid/ChIP_IgG-100221_S18.sort.rmdup.rmblackls.rmchr.bam"

#################### MAIN CODE ###########################

# full path of bams to subsample
tosubsample=${outPath}"/subsampling/toSubsample.txt"

outbampath=${outPath}"/subsampling/bamfiles"
outpeak=${outPath}"/subsampling/peaks"
summaryPath=${outPath}"/subsampling/summary"

if [ ! -e ${outbampath} ]; then
        mkdir -p ${outbampath}
fi
if [ ! -e ${outpeak} ]; then
        mkdir -p ${outpeak}
fi
if [ ! -e ${summaryPath} ]; then
        mkdir -p ${summaryPath}
fi


#inbam="/home/jmendietaes/data/2021/chip/allProcessed/bamfiles/valid/Mye_Smarcb1_31032021_S23.sort.rmdup.rmblackls.rmchr.bam"

for inbam in `cat ${tosubsample}`; do

    inbam_name=`basename ${inbam}`
    

    if [[ $inbam_name == Mye* && $inbam_name != $(basename ${myeControl}) ]] ; then
        controlbam=${myeControl}
    elif [[ $inbam_name == DM* && $inbam_name != $(basename ${dmControl}) ]] ; then
        controlbam=${dmControl}
    elif [[ $inbam_name == ChIP* && $inbam_name != $(basename ${chipControl}) ]] ; then
        controlbam=${chipControl}
    else
        controlbam=""
    fi


    ###############################
    # peak calling of whole file
    ###############################
    # Create output dir

    if [ ! -e ${outpeak}/logs ]; then
        mkdir -p ${outpeak}/logs
    fi

    if [ ! -e ${summaryFile} ] ; then
        touch ${summaryFile}
    fi


    ip=`echo $inbam_name | sed 's/\.sort.*//g'`
    summaryFile="${summaryPath}/${ip}.summary"

    echo -e "Starting Peak calling of full bam file -------------------------------------- \n"

    # check content of first line of step control file
    if grep -q "NUMBER OF NARROW PEAKS$(printf '\t')1" ${summaryFile}; then
        echo "Peak calling of full bam already done before ----------------------- \n"
    else
        total_reads=$(samtools view -c ${inbam})

        # narrow peaks
        peaktype='narrowPeak'
        macs2 callpeak \
                -t ${inbam} \
                -c ${controlbam} \
                -f BAMPE \
                -g $species \
                -n $ip \
                --keep-dup all \
                --outdir ${outpeak}/ 2> ${outpeak}/logs/${ip}_macs2.log

        mv ${outpeak}/${ip}_peaks.xls ${outpeak}/${ip}_peaks_${peaktype}.xls

        npeaks=$(cat ${outpeak}/${ip}_peaks.${peaktype} | wc -l)
        reads_in_peaks=$(bedtools sort -i ${outpeak}/${ip}_peaks.${peaktype} \
            | bedtools merge -i stdin | bedtools intersect -u -nonamecheck \
            -a ${inbam} -b stdin -ubam | samtools view -c)
        FRiP=$(awk "BEGIN {print "${reads_in_peaks}"/"${total_reads}"}")
        # report
        echo -e "NUMBER OF NARROW PEAKS\t1\t${npeaks}" >> ${summaryFile}
        echo -e "total_reads\treads_in_peaks\tFRIP" >> ${summaryFile}
        echo -e "${total_reads}\t${reads_in_peaks}\t${FRiP}\n" >> ${summaryFile}

        echo -e "Peak calling of full bam - Done ------------------------- \n"
    fi






    ###################################
    # subsample and new peak callings
    ###################################
    for fractionOfReads in `seq 0.1 0.2 1`; do
        echo -e "Starting Peak calling of ${fractionOfReads} file ------------------------- \n"

        if grep -q "NUMBER OF NARROW PEAKS$(printf '\t')${fractionOfReads}" ${summaryFile}; then
            echo "Peak calling of ${fractionOfReads} already done before ----------------------- \n"
        else

            subSbam=`echo $inbam_name | sed "s/rmblackls\.rmchr/rmblackls\.rmchr_${fractionOfReads}/g"`
            subSbam="${outbampath}/${subSbam}"

            # subsample and create bam index
            sambamba view -h -t $SLURM_CPUS_PER_TASK -s $fractionOfReads -f bam --subsampling-seed=12345 $inbam -o $subSbam


            ip2="${ip}_${fractionOfReads}"
            total_reads=$(samtools view -c ${subSbam})

            # narrow peaks

            peaktype='narrowPeak'
            if [[ $controlbam == "" ]] ; then
                echo -e "No control file used"
                macs2 callpeak \
                    -t ${subSbam} \
                    -f BAMPE \
                    -g $species \
                    -n $ip2 \
                    --keep-dup all \
                    --outdir ${outpeak}/ 2> ${outpeak}/logs/${ip2}_macs2.log
            else
                macs2 callpeak \
                        -t ${subSbam} \
                        -c ${controlbam} \
                        -f BAMPE \
                        -g $species \
                        -n $ip2 \
                        --keep-dup all \
                        --outdir ${outpeak}/ 2> ${outpeak}/logs/${ip2}_macs2.log
            fi

            mv ${outpeak}/${ip2}_peaks.xls ${outpeak}/${ip2}_peaks_${peaktype}.xls

            npeaks=$(cat ${outpeak}/${ip2}_peaks.${peaktype} | wc -l)
            reads_in_peaks=$(bedtools sort -i ${outpeak}/${ip2}_peaks.${peaktype} \
                | bedtools merge -i stdin | bedtools intersect -u -nonamecheck \
                -a ${subSbam} -b stdin -ubam | samtools view -c)
            FRiP=$(awk "BEGIN {print "${reads_in_peaks}"/"${total_reads}"}")
            # report
            echo -e "NUMBER OF NARROW PEAKS\t${fractionOfReads}\t${npeaks}" >> ${summaryFile}
            echo -e "total_reads\treads_in_peaks\tFRIP" >> ${summaryFile}
            echo -e "${total_reads}\t${reads_in_peaks}\t${FRiP}\n" >> ${summaryFile}

            echo -e "Peak calling of ${fractionOfReads} file - Done ------------------------- \n"
            rm ${subSbam}
        fi
    done
done