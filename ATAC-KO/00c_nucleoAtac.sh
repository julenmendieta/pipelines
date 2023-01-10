# HOW TO RUN ME
bash /home/julen/programas/PhD/ATAC-KO/00c_nucleoAtac.sh

# This scripts runs TOBIAS comparing datasets with a specific
# string with others. It also discerns comparison by batches
# specified in the name after the first _
# cell-condition_batch-[extra]_...

# Path to nucleoAtac
nucleoatac=/home/julen/miniconda3/envs/python2/bin/nucleoatac
# path to genome fasta file
fastap=/home/julen/genomes/mm10_reordered/mm10.reordered.fa
# Number of cores to use
nCores=16
# base bams path
bambase=/scratch/julen/ATAC/allData/02_firstATAC/outdata/bams/merged
# base peaks path
peakbase=/scratch/julen/ATAC/allData/02_firstATAC/peaks
# define peak type
peaktype=narrowPeak
# path to main output folder (ending with '/')
outbase=/scratch/julen/ATAC/allData/02_firstATAC/nucleoAtac

##===============================================================================

# function to check if the given first file doesnt exist or is older than 
# the second input file
fileNotExistOrOlder () {
    # check if the file exists of it was created with a previous bam version 
    analyse="no"
    if [ ! -e $1 ]; then
        analyse="yes"
    # only proceed if the output file is older than the bam file
    # in this way if we resequenced and kept the name the analysis 
    # will be repeated
    else
        for tfile in $2; do
            if [[ $1 -ot ${tfile} ]] ; then
                analyse="yes"
                echo $1" older than "${tfile}
            fi
        done
    fi
}

##===============================================================================



nucleoatac run --bed <bedfile> --bam <bamfile> --fasta ${fastap} \
            --out <output_basename> --cores ${nCores}


bedfile=${peakbase}/LSK-NTC_ATAC2-merged_peaks.broadPeak
bamfile=${bambase}/LSK-NTC_ATAC2-merged.sort.rmdup.rmblackls.rmchr.Tn5.bam

$nucleoatac run --bed ${bedfile} --bam ${bamfile} --fasta ${fastap} \
            --out ${outbase}/LSK-NTC_ATAC2-merged --cores ${nCores}

# Convert bedgraphs to bigwigs
gzip -d nucleoAtac.nucleoatac_signal.bedgraph.gz
bedGraphToBigWig nucleoAtac.nucleoatac_signal.bedgraph /home/julen/genomes/mm10_reordered/mm10.reordered.sizes.short nucleoatac_signal.bedgraph.bw



chip=LSK-Smarcd2_ATAC2-merged
bedfile=${peakbase}/${chip}_peaks.broadPeak
bamfile=${bambase}/${chip}.sort.rmdup.rmblackls.rmchr.Tn5.bam

$nucleoatac run --bed ${bedfile} --bam ${bamfile} --fasta ${fastap} \
            --out ${outbase}/${chip} --cores ${nCores}

