#!/bin/bash
# -*- ENCODING: UTF-8 -*-

# HOW TO RUN ME
# Change the variables to change and
# bash /home/jmendietaes/programas/PhD/general/renameFiles.sh

# OBJECTIVE OF SCRIPT
# To read a file with two columns, and convert the file names 
#in given path that contain the string in column 1 by 
#substituting the string by the one from column 2

########## TO CHANGE ###########
# Path to the table with 2 columns (string\tsubstitute)
# DO NOT USE home relative paths (~/)
inTable="/home/jmendietaes/data/2021/singleCell/allProcessed/pruebas/convertNames.txt"
# path in which we will do the search and name changing
focusPath="/home/jmendietaes/data/2021/singleCell/allProcessed/cloupes"
# If we want to do a test printing the change but not doing anything
# posible answers are lowercase "yes" (ony print) or "no" (print and change)
onlyTest="no"


########## START #############
# store table content in arrays (strings and substitutes)
IFS=$'\n' read -r -d '' -a strings < <( awk '{print $1}' ${inTable} && printf '\0' )
IFS=$'\n' read -r -d '' -a substitutes < <( awk '{print $2}' ${inTable} && printf '\0' )

# make sure they have same length
if [ ! "${#strings[@]}" -eq "${#substitutes[@]}" ]; then
    echo "Different number of elements in each column"
    exit 1
fi

for ni in "${!strings[@]}"; do
    string_=${strings[${ni}]}
    substitute_=${substitutes[${ni}]}
    file=$(find ${focusPath} -maxdepth 1 -name *${string_}* | grep .)

    # if we found a match
    if [ ! -z "$file" ]; then 
        for fi in $file; do 
            basefi=$(basename ${fi})
            newName=$(echo ${basefi} | sed "s/${string_}/${substitute_}/g")
            echo "$basefi   ---    $newName"
            if [ ${onlyTest} == "no" ]; then 
                mv ${focusPath}/${basefi} ${focusPath}/${newName}
            fi
        done
    fi
done
