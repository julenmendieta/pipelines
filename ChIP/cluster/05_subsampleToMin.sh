#!/bin/bash
# -*- ENCODING: UTF-8 -*-

##===============================================================================
## SLURM VARIABLES
#SBATCH --job-name=SubsToMin
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=05:00:00
#SBATCH -p short
#SBATCH -o /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.out  
#SBATCH -e /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.err 

# HOW TO RUN ME
#sbatch /home/jmendietaes/programas/PhD/ChIP/cluster/05_subsampleToMin.sh \
#/home/jmendietaes/data/2021/chip/allProcessed \

# PURPOSE
# to subsample Mye or DM bam files to the one that has less reads to then
# store it in the subsampled folder

basePath=$1
#basePath="/home/jmendietaes/data/2021/chip/allProcessed"

# cells to be counted at the time to subsample (separated by \|)
subsampleCells="DM_\|Mye_"

bamsPath="${basePath}/bamfiles/valid"
subOut="${bamsPath}/subsampled"
bamCounts="${bamsPath}/bamCounts.txt"

# modules
module load Sambamba/0.7.0

# Create output dir
if [ ! -e ${subOut} ]; then
    mkdir -p ${subOut}
fi

# function to join elements from array
function join_by { local IFS="$1"; shift; echo "$*"; }

# Make sure bamCounts.txt is updated
# this will overwrite it (very slow, better if you update what you need)
#echo '' > ${bamCounts}
cd ${bamsPath}
checkBams=$(ls *bam | grep "${subsampleCells}") # | grep "Brd9\|CTCF\|EHMT1\|RBBP4")
for i in ${checkBams}; do 
    fileSize=$(du -k ${i} | cut -f1); 
    fileSize_prev=$({ grep $i ${bamCounts} | cut -f 1 || :; });
    # if first time we check it
    if [[ ${fileSize_prev} == "" ]] ; then
        echo $i;
        fileLen=$(samtools idxstats ${i} | cut -f 3 | awk '{sum+=$1;} END{print sum;}');
        echo -e "${fileSize}\t${i}\t${fileLen}" >> bamCounts.txt ; 
    elif [[ ${fileSize} != ${fileSize_prev} ]]; then
        echo $i;
        # store file without this entry
        grep -v $i ${bamCounts} > bamCounts.txt_ ;
        mv bamCounts.txt_ bamCounts.txt ;
        fileLen=$(samtools idxstats ${i} | cut -f 3 | awk '{sum+=$1;} END{print sum;}');
        echo -e "${fileSize}\t${i}\t${fileLen}" >> bamCounts.txt ; 
    fi;
done ; for i in mergedReplicates/*bam; do 
    fileSize=$(du -k ${i} | cut -f1); 
    fileSize_prev=$({ grep $i ${bamCounts} | cut -f 1 || :; });
    # if first time we check it
    if [[ ${fileSize_prev} == "" ]] ; then
        echo $i;
        fileLen=$(samtools idxstats ${i} | cut -f 3 | awk '{sum+=$1;} END{print sum;}');
        echo -e "${fileSize}\t${i}\t${fileLen}" >> bamCounts.txt ; 
    elif [[ ${fileSize} != ${fileSize_prev} ]]; then
        echo $i;
        # store file without this entry
        grep -v $i ${bamCounts} > bamCounts.txt_ ;
        mv bamCounts.txt_ bamCounts.txt ;
        fileLen=$(samtools idxstats ${i} | cut -f 3 | awk '{sum+=$1;} END{print sum;}');
        echo -e "${fileSize}\t${i}\t${fileLen}" >> bamCounts.txt ; 
    fi;
done

chips=$(for fi in `cut -f 2 ${bamCounts} | tail -n +2`; do 
        filename=$(basename ${fi}); 
        mapLib=(${filename//\./ }); 
        mapLib=${mapLib[0]};
        mapLib=(${mapLib//_/ }); 
        mapLib=${mapLib[1]}; 
        mapLib=(${mapLib//-/ }) ; 
        echo ${mapLib[0]};
    done | sort| uniq)

# To show files and number of reads
#for chip in $chips; do echo $chip; grep $chip ${bamCounts} | cut -f 2,3; echo ; done

# for now we only subsample DM and Mye
for chip in $chips; do
    echo $chip
    chipFiles=$(grep "${chip}_\|${chip}-\|${chip}\." ${bamCounts}| cut -f 2 | \
                    grep -v mergedReplicates | grep "${subsampleCells}")

    if [[ $(echo $chipFiles | wc -w) == 2 ]]; then
        # get minimum number of reads
        minRead=$(for chi in $chipFiles; do grep $chi ${bamCounts} | cut -f 3; done | sort -n | head -n 1)

        # define new names per file 
        for file1 in ${chipFiles}; do
            mapLib=(${file1//_/ }); 
            # we only had one _
            if [[ ${#mapLib[@]} == 2 ]]; then
                mapLib=(${file1//\./ }); 
                mapLib[0]=${mapLib[0]}"-sub${minRead}"
                subsName=$(join_by . ${mapLib[@]})
            # we had more than one
            else
                mapLib[1]=${mapLib[1]}"-sub${minRead}"
                subsName=$(join_by _ ${mapLib[@]})
            fi

            # get proportion and subsample
            nReads=$(grep $file1 ${bamCounts} | cut -f 3)
            if [[ ${nReads} == ${minRead} ]]; then
                # if file doesnt exist
                if [ ! -e ${subOut}/${subsName} ]; then
                    cp ${bamsPath}/${file1} ${subOut}/${subsName}
                    cp ${bamsPath}/${file1}.bai ${subOut}/${subsName}.bai
                fi
            else
                # if file doesnt exist
                if [ ! -e ${subOut}/${subsName} ]; then
                    fractionOfReads=$(echo "print(${minRead}/${nReads})" | python3)
                    sambamba view -h -t $SLURM_CPUS_PER_TASK -s $fractionOfReads \
                            -f bam --subsampling-seed=12345 ${bamsPath}/${file1} \
                            -o ${subOut}/${subsName}
                fi
                
            fi


        done
    fi
done



# # To copy only subsampled and not subsampled files that are not common
# chips=$(for fi in `cat ${bamCounts} | grep "${subsampleCells}" | \
#         grep -v mergedReplicates/ | cut -f 2  | tail -n +2`; do 
#         filename=$(basename ${fi}); 
#         mapLib=(${filename//\./ }); 
#         mapLib=${mapLib[0]};
#         mapLib=(${mapLib//_/ }); 
#         mapLib=${mapLib[1]}; 
#         mapLib=(${mapLib//-/ }) ; 
#         echo ${mapLib[0]};
#     done | sort| uniq)

# chips=$(for chip in $chips; do
#     chipFiles=$(grep "${chip}_\|${chip}-" ${bamCounts}| cut -f 2 | \
#                     grep -v mergedReplicates | grep "${subsampleCells}")

#     if [[ $(echo $chipFiles | wc -w) == 2 ]]; then
#         echo "$chip"
#     fi
# done)
# cd /home/jmendietaes/data/2021/chip/allProcessed/bamfiles/valid/08_projectRestart_subsampled
# for chip in ${chips}; do 
#     ln -s /home/jmendietaes/data/2021/chip/allProcessed/bamfiles/valid/subsampled/DM_${chip}* .
#     ln -s /home/jmendietaes/data/2021/chip/allProcessed/bamfiles/valid/subsampled/Mye_${chip}* .
# done

# vgrep=$(for chip in $chips; do echo "${chip}\|"; done | tr '\n' ' ')
# vgrep=$(echo $vgrep | sed 's/ //g')
# vgrep=${vgrep::-2}
# for i in `ls /home/jmendietaes/data/2021/chip/allProcessed/bamfiles/valid/DM_* | grep -v "${vgrep}"`; do
#     ln -s ${i} .
# done
# for i in `ls /home/jmendietaes/data/2021/chip/allProcessed/bamfiles/valid/Mye_* | grep -v "${vgrep}"`; do
#     ln -s ${i} .
# done