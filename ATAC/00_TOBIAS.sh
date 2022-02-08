

# Path to TOBIAS
TOBIAS=/home/julen/miniconda3/envs/TOBIAS/bin/TOBIAS
# base bams path
bambase=/scratch/julen/ATAC/bams/merged
# base peaks path
peakbase=/scratch/julen/ATAC/allData/02_firstATAC/peaks
# base output dir
outdir=/scratch/julen/ATAC/allData/02_firstATAC/TOBIAS
# Path to motifs to bhe check
motifsCheck=~/programas/HOMER/data/knownTFs/vertebrates/known.motifs
# Path to reference genome
refGenome=/home/julen/genomes/mm10_reordered/mm10.reordered.fa
# Cell type of interest
cell=LSK
# define peak type
peaktype=narrowPeak
# path for the location of the pipeline scripts
scriptsPath="/home/julen/programas/PhD"
# Number of CPU to use
nCPU=16



# Convert homer motifs to JASPAR  matrix
file="/scratch/julen/ATAC/allData/02_firstATAC/TOBIAS/20220204133836_JASPAR2022_combined_matrices_17766_pfm.txt"
fileout=/scratch/julen/ATAC/allData/02_firstATAC/TOBIAS/motifClustering/topHomer/allTogether_matrix.motifs
motifs = read_homer(file, skip = 0)
motifsPFM = convert_motifs(motifs,class = "motifStack-pfm")
write_matrix(motifs, fileout, positions = "columns", rownames = FALSE,
        'jaspar', sep = "", headers = TRUE, overwrite = TRUE,
        append = FALSE)


if [ ! -e ${outdir}/${cell} ]; then
	mkdir -p ${outdir}/${cell}
fi

# get NTC and all KO
allbams=$(find -L ${bambase}/${cell}*bam)
NTC_bam=$(echo $allbams| tr ' ' '\n' | grep "${cell}-NTC_") 
ko_bams=$(echo $allbams| tr ' ' '\n' | grep -v "${cell}-NTC_") 

#########################
# CONSENSUS PEAKS: Get consensus peaks for same cell type
#########################

# Lets create a consensus peak coordinates file for all the conditions of same cell
allPeaks=$(find -L ${peakbase}/${cell}*${peaktype})
if [ ${peaktype} == "narrowPeak" ]; then
    mergecols=`seq 2 10 | tr '\n' ','`
    expandparam='--is_narrow_peak'
elif [ ${peaktype} == "broadPeak" ]; then
    mergecols=`seq 2 9 | tr '\n' ','`
    expandparam=''
fi
fileLabels=$(for f in $allPeaks; do echo ${f##*/} |sed "s/_peaks.${peaktype}//g"; done | tr '\n' ',')

sort -T '.' -k1,1 -k2,2n ${allPeaks} \
                | mergeBed -c $mergecols -o collapse > ${outdir}/${cell}/${cell}_consenus.txt
python ${scriptsPath}/ChIP/cluster/02_NR_macs2_merged_expand.py ${outdir}/${cell}/${cell}_consenus.txt \
                    ${fileLabels} \
                    ${outdir}/${cell}/${cell}_consenus.boolean.txt \
                    $expandparam

consensusPeakBed=${outdir}/${cell}/${cell}.bed
awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print $1, $2, $3, $4, "0", "+" }' \
    ${outdir}/${cell}/${cell}_consenus.boolean.txt > ${consensusPeakBed}
rm ${outdir}/${cell}/${cell}_consenus.txt
rm ${outdir}/${cell}/${cell}_consenus.boolean.txt



#########################
# NTC analysis: Run first steps in NTC
#########################

prefix_ntc=$(basename ${NTC_bam} | sed 's/.sort.rmdup.*.bam//g')
outNTC=${outdir}/${cell}/${prefix_ntc}
bamfile=${bambase}/${prefix_ntc}.sort.rmdup.rmblackls.rmchr.Tn5.bam

$TOBIAS ATACorrect --bam ${bamfile} \
        --genome ${refGenome} \
        --peaks ${consensusPeakBed} --read_shift 0 0 --prefix ${prefix_ntc} \
        --outdir ${outNTC} --cores ${nCPU}

#Calculate footprint scores per condition*/
correctedbw=${outNTC}/${prefix_ntc}_corrected.bw
$TOBIAS ScoreBigwig --signal ${correctedbw} --regions ${consensusPeakBed} \
                    --cores ${nCPU} --output ${outNTC}/${prefix_ntc}_footprints.bw &> \
                    ${outNTC}/${prefix_ntc}_footprints.log
    

#########################
# KO analysis: Run steps in KO
#########################

# Run steps and KO samples and then get differential sites
for ko_bam in ${ko_bams}; do
    prefix_ko=$(basename ${ko_bam} | sed 's/.sort.rmdup.*.bam//g')
    outKo=${outdir}/${cell}/${prefix_ko}
    bamfile=${bambase}/${prefix_ko}.sort.rmdup.rmblackls.rmchr.Tn5.bam

    $TOBIAS ATACorrect --bam ${bamfile} \
            --genome ${refGenome} \
            --peaks ${consensusPeakBed} --read_shift 0 0 --prefix ${prefix_ko} \
            --outdir ${outKo} --cores ${nCPU}
    #TOBIAS ATACorrect ${atacorrect} -b $allmerged -g $fasta -p $allbed --cores 99 ${atacorrect} --blacklist $blacklist --prefix $condition &> ${condition}_atacorrect.log
    #Calculate footprint scores per condition*/
    correctedbw=${outKo}/${prefix_ko}_corrected.bw
    $TOBIAS ScoreBigwig --signal ${correctedbw} --regions ${consensusPeakBed} \
                    --cores ${nCPU} --output ${outKo}/${prefix_ko}_footprints.bw &> \
                    ${outKo}/${prefix_ko}_footprints.log
    #TOBIAS ScoreBigwig --signal ${correctedbw} ${footprinting} --regions $allmerged ${footprinting} --cores 99 --output ${condition}_footprints.bw &> ${condition}_footprinting.log

    #Estimate bound sites from scored */
    # here we also can do differential
    
    outprefix1=(${prefix_ntc//_/ }); 
    outprefix1=${outprefix1[0]};
    outprefix2=(${prefix_ko//_/ }); 
    outprefix2=${outprefix2[0]}; 
    outprefix=${outprefix1}_${outprefix2}

    # It works as it should when you use raw PFM motif files
    motifsCheck=/scratch/julen/ATAC/allData/02_firstATAC/TOBIAS/motifClustering/topHomer/allTogether.motifs

    $TOBIAS BINDetect --motifs ${motifsCheck} \
                    --signals ${outNTC}/${prefix_ntc}_footprints.bw ${outKo}/${prefix_ko}_footprints.bw \
                    --genome ${refGenome} --peaks ${consensusPeakBed} \
                    --cond-names ${prefix_ntc} ${prefix_ko} --cores ${nCPU} \
                    --prefix bindetect-${outprefix} \
                    --debug \
                    --outdir ${outKo}/bindetect_output &> ${outKo}/${outprefix}_bindetect.log
    #TOBIAS BINDetect --motifs ${motifsCheck} --signals $footprints ${bindetect} --genome $fasta --peaks  $annotated_headerbed --peak_header $annotated_header --cores 99 --cond_names $condition --outdir TFBS &> bindetect.log
    #Join subset of bound TFs*/
    cat $bed | cut -f1-4 | sort -k1,1 -k2,2n > all_${condition}_bound.bed; igvtools index all_${condition}_bound.bed 2>&1;


## Appart from all
$TOBIAS ClusterMotifs --motifs ${motifsCheck} --threshold 0.3 --type png --outdir /scratch/julen/ATAC/allData/02_firstATAC/TOBIAS/motifClustering