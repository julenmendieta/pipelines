

inpath=/home/jmendietaes/data/dataSubmission/2022-ainhoaPaper/ChIP/bamfiles
filePattern="bam"
geoIN=julenm_OZS9ouGE/geo_submission_sept13/ChIP
LOGIN="geoftp"
#PW=$1
PW=rebUzyi1



#########################################
# exit when any command fails
#set -e

logfile="${inpath}/Files_Transfer.tsv" # in this file the names of the files to transfer are stored
touch $logfile
# Get files where upload is completed
filesCompleted=$(lftp -u $LOGIN,$PW ftp-private.ncbi.nlm.nih.gov/uploads/${geoIN} -e "ls; exit")
filesCompleted=$(echo $filesCompleted | tr ' ' '\n' | grep ${filePattern})
echo $filesCompleted

# upload other files
files=$(find $inpath/*${filePattern} -maxdepth 0 -mindepth 0)
#line="$GFS/PROJECTS/JAKSTAT/Results_pipelines/aligned_bam/RNA_a8_ko_t_1.Aligned.out_sorted.bam RNA_T8_STAT6KO_H_H_SKAF001.bam"
#while read line; do
for fi in ${files}; do
  echo "---------------------------------"
  #echo "$fi"
  #originalPath=$(echo $line | sed "s/\s.*//")
  targetName=$(basename $fi)
  #targetName=$(echo $line | sed "s/.*\s//")
  echo "-- Files:"
  #echo $originalPath
  #echo $originalFile
  echo $targetName
 
  echo "-- Transfer:"
  if echo $filesCompleted | grep -q "$targetName"; then
    echo "$targetName transferred"
  else
    echo "Transfering $targetName"
    lftp -u $LOGIN,$PW ftp-private.ncbi.nlm.nih.gov/uploads/${geoIN} -e "put $targetName; exit"
    #lftp -u $LOGIN,$PW ftp-private.ncbi.nlm.nih.gov/uploads/${geoIN} -e "mv $originalFile $targetName; exit"

  fi
done >> $logfile