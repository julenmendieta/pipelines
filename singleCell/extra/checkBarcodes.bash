for i in *_*; do
    nfiles=$(ls ${i} | wc -w)
    if [[ ${nfiles} -ge 2 ]]; then 
        echo
        echo $i
        for dataT in mRNA gRNA; do
            echo
            echo ${dataT}
            for run in `ls -d $i/*`; do
                if [ -e ${run}/${dataT} ]; then
                    basename $run
                    zcat ${run}/${dataT}/*_R1*gz | head -n 1
                fi
            done
        done
    fi
done
