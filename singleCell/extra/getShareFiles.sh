# Script to get and copy cloupe and stats files of CellRanger output 
# (Will rename them to contain sample ID)
# cloupe files
echo "Normal files"
for i in /home/jmendietaes/data/2021/singleCell/allProcessed/Data/*; do
    if [ -e "${i}/outs" ]; then
        basei=$(basename $i)
        if [[ $i == *onlyRNA ]] ; then
            echo "This file shouldnt be here"
            echo $i
            #if [ ! -f "cloupes/onlyRNA/$basei.cloupe" ] ; then
            #    cp ${i}/outs/cloupe.cloupe cloupes/onlyRNA/$basei.cloupe;
            #fi
        else
            if [ ! -f "cloupes/$basei.cloupe" ] ; then
                echo $i
                cp ${i}/outs/cloupe.cloupe cloupes/$basei.cloupe;
            fi
        fi
    fi
done

#echo "onlyRNA"
#for i in /home/jmendietaes/data/2021/singleCell/allProcessed/Data/onlyRNA/*; do
#    if [ -e "${i}/outs" ]; then
#        basei=$(basename $i)
#        if [[ $i == *onlyRNA ]] ; then
#            if [ ! -f "cloupes/onlyRNA/$basei.cloupe" ] ; then
#                echo $i
#                cp ${i}/outs/cloupe.cloupe cloupes/onlyRNA/$basei.cloupe;
#            fi
#        else
#            echo "This file shouldnt be here"
#            echo $i
#            #if [ ! -f "cloupes/$basei.cloupe" ] ; then
#            #    cp ${i}/outs/cloupe.cloupe cloupes/$basei.cloupe;
#            #fi
#        fi
#    fi
#done

################
# stats files
echo "Normal files"
for i in /home/jmendietaes/data/2021/singleCell/allProcessed/Data/*; do
    if [ -e "${i}/outs" ]; then
        basei=$(basename $i)
        if [[ $i != *onlyRNA ]] ; then
            if [ ! -f "stats/${basei}_web_summary.html" ] ; then
                echo $i
                cp ${i}/outs/web_summary.html stats/${basei}_web_summary.html;
            fi
        fi
    fi
done

#echo "onlyRNA"
#for i in /home/jmendietaes/data/2021/singleCell/allProcessed/Data/onlyRNA/*; do
#    if [ -e "${i}/outs" ]; then
#        basei=$(basename $i)
#        if [[ $i == *onlyRNA ]] ; then
#            if [ ! -f "stats/onlyRNA/${basei}_web_summary.html" ] ; then
#                echo $i
#                cp ${i}/outs/web_summary.html stats/onlyRNA/${basei}_web_summary.html;
#            fi
#        fi
#    fi
#done
