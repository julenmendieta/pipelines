# TO BE EXECUTED LINE BY LINE

# https://metascape.org/gp/index.html#/menu/msbio

# Start and stop docker
#systemctl start docker
#systemctl stop docker
#systemctl stop docker.socket

# Path to MSBio container
msPath="/home/julen/programas/MSBio"
# Set to "-u" if input format follows the single-gene-list standard (else to "")
geneListStyle=""
# Proyect name (input and output folders are inside ${msPath}/data/${pName})
pName=KOtogether_overlapOnly
# Place your input txt files in ${msPath}/data/${pName}/input

nCPU=8
# For -S and -T you can use either taxonomy ID or common names. 
# The supported IDs are: 9606, 10090, 10116, 4932, 5833, 6239, 7227, 7955, 
#   3702, and 4896. 
# The supported names are: human, mouse, rat, yeast, malaria, "c. elegans", 
#   fly, zebrafish, arabidopsis, or "s. pombe".
# If the gene list is not for human (default), use -S SOURCE_TAX_ID
taxId=10090
# If the target organism is not for human (default), use -T TARGET_TAX_ID
targetId=9606

cd ${msPath}
# Start MSBio
sudo ${msPath}/bin/up.sh

# If we want to run all input files
inputDir="${msPath}/data/${pName}/input"
optionJson="/data/${pName}/option.json"
#optionJson="/data/option.json"
files=$(ls $inputDir)
for input_list_name in $files; do
    echo $input_list_name
    id1=$(echo $input_list_name | sed 's/\.txt//g')

    # The input and output folders must be subfolders of ${msPath}/data
    # This folder is mounted in the container as /data
    output_folder="/data/${pName}/${id1}"
    input_list_file="/data/${pName}/input/${id1}.txt"

    # Analyse gene list
    # You must specify -u if your input format follows the single-gene-list standard
    sudo ${msPath}/bin/ms.sh ${singleGeneList} -o ${output_folder} \
                            -S ${taxId} -T ${targetId} ${input_list_file} \
                            --option ${optionJson} -c ${nCPU} ${geneListStyle}

done

# End MSBio
sudo ${msPath}/bin/down.sh


# # Version for a specific file
# # The input and output folders must be subfolders of ${msPath}/data
# # This folder is mounted in the container as /data
# output_folder="/data/${pName}/Mye_100000bp"
# input_list_file="/data/${pName}/input/Mye_clusterGenes_100000bp_plusBg.txt"

# # Start MSBio
# sudo ${msPath}/bin/up.sh

# # Analyse gene list
# # You must specify -u if your input format follows the single-gene-list standard
# sudo ${msPath}/bin/ms.sh ${singleGeneList} -o ${output_folder} \
#                         -S ${taxId} -T ${targetId} ${input_list_file}

# # End MSBio
# sudo ${msPath}/bin/down.sh



## NOTES
# Although Metascape prioritizes results for terms that are commonly 
#   shared across gene lists, users can check “Pick selective GO 
#   clusters” to adopt a different prioritization algorithm to 
#   preferably identify terms that are selective across gene lists 
#   instead. (l_go_selective ?)

# How to copy zip links
#files=$(ls -d ~/programas/MSBio/data/perCell_expreUniq/*_*)
#for fi in $files; do fname=$(basename $fi); ln -s ${fi}/all.zip ${fname}.zip; done
# how to unzip
#for zipf in *zip; do bname=$(echo $zipf | sed "s/\.zip//g"); unzip ${zipf} -d ${bname}; done
