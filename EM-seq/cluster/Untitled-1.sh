module load SAMtools/1.12-GCC-10.2.0
module load BEDTools/2.27.1-foss-2018b
module load zlib/1.2.11-GCCcore-11.2.0
export PATH="/home/jmendietaes/programas/miniconda3/envs/DNAme/bin:$PATH"

nCPU=8
PROJECT_DIR="/home/jmendietaes/data/2021/DNAme/sequencedData/NextSeq2000.RUN156.20230306"
filename=DM_DNAme1_1_S22
RAW_FASTQ_DIR=$PROJECT_DIR"/demux_fastq"
REFERENCE_DIR="/home/jmendietaes/referenceGenomes/mm10_reordered/mm10.reordered"
read1_path="${RAW_FASTQ_DIR}/${filename}_R1_001.fastq.gz"
read2_path="${RAW_FASTQ_DIR}/${filename}_R2_001.fastq.gz"
bamFile="${EDITED_DIR}/BAM/${filename}.sort.markdup.bam"
EDITED_DIR=$PROJECT_DIR"/pipelineOut"
basePath="/home/jmendietaes/data/2021/DNAme/allProcessed"
summaryFile="${basePath}/QC/summary_${filename}.txt"

# Index reference genome
# biscuit index my_reference.fa

# Align sequencing reads to the reference
# Gzipped FASTQ files can also be used
biscuit align -@ ${nCPU} ${REFERENCE_DIR}.fa ${read1_path} ${read2_path} |
    dupsifter ${REFERENCE_DIR}.fa | samtools sort -@ ${nCPU} -o my_output.bam -O BAM -
samtools index -@ ${nCPU} my_output.bam

# Run QC
biscuit qc /path/to/my_reference.fa input.bam sample_name
${subScripts}/QC.sh [--vcf input.vcf] \
    [--outdir output_dir] \
    [--no-cov-qc] \
    /home/jmendietaes/referenceGenomes/mm10_reordered/Biscuit/assets/mm10_biscuit_qc_assets.zip \
    /path/to/my_reference.fa \
    sample_name \
    input.bam



# Filter Reads by Bisulfite Conversion
# The -p flag outputs the counts in a table, instead of as a tag in the BAM file
biscuit bsconv /path/to/my_reference.fa input.bam [out.bam]

# Validate Bisulfite Conversion Label
biscuit bsstrand [-c] /path/to/my_reference.fa input.bam [out.bam]
# f: OT/CTOT (BSW) strand
# r: OB/CTOB (BSC) strand
# c: conflicting strand information
# u: unintelligible strand source (unknown)


# Create a pileup VCF of DNA methylation and genetic information
# Also compresses and indexes the VCF
# allows the user to compute cytosine retention and callable SNP mutations.
# BISCUIT has the ability to put mutation calls and DNA methylation 
# measurements from multiple samples next to each other in the output VCF 
# by providing biscuit pileup with more than one input BAM
biscuit pileup -@ ${nCPU} -o my_pileup.vcf ${REFERENCE_DIR}.fa my_output.bam
bgzip -@ ${nCPU} my_pileup.vcf
tabix -p vcf my_pileup.vcf.gz

# Extract DNA methylation into BED format
# Also compresses and indexes the BED
biscuit vcf2bed -t cg my_pileup.vcf.gz > my_methylation_data.bed
bgzip my_methylation_data.bed
tabix -p bed my_methylation_data.bed.gz

#  Merge neighboring C and G in CpG context
biscuit mergecg /path/to/my_reference.fa my_methylation_data.bed

# Get file for Bismar COV (Methylation coverage track)
biscuit vcf2bed -c my_pileup.vcf.gz > my_beta_m_u.bed
awk -v OFS='\t' '{ print $1, $2+1, $3, $4, $5, $6 }' my_beta_m_u.bed > my_beta_m_u.cov
# This file will keep methylation across the strands separate, as is normally 
# done in biscuit vcf2bed. To merge methylation across strands, run:
#biscuit vcf2bed my_pileup.vcf.gz | \
#biscuit mergecg -c /path/to/my_reference.fa - | \
#awk -v OFS='\t' '{ print $1, $2+1, $3-1, $4, $5, $6 }' > my_merged_beta_m_u.cov





This basic order of commands will produce all the necessary files needed to 
read data into R using the R/Bioconductor companion package, biscuiteer.
https://www.bioconductor.org/packages/release/bioc/html/biscuiteer.html


# BISCUIT has the ability to put mutation calls and DNA methylation 
# measurements from multiple samples next to each other in the output VCF 
# by providing biscuit pileup with more than one input BAM
biscuit pileup -@ ${nCPU} -o my_combined_pileup.vcf ${REFERENCE_DIR}.fa my_outpu1t.bam my_output2.bam
bgzip -@ ${nCPU} my_combined_pileup.vcf
tabix -p vcf my_combined_pileup.vcf.gz

# BISCUIT can call somatic mutations by providing a tumor and matched normal
#  BAM to pileup. To run in somatic mode, run biscuit pileup with the -S flag:

# biscuit pileup -@ NTHREADS -S -o somatic_mode.vcf /path/to/my_reference.fa \
#     -T tumor.bam -I normal.bam
# bgzip -@ NTHREADS somatic_mode.vcf
# tabix -p vcf somatic_mode.vcf.gz


Pileup VCF files were generated using biscuit pileup with default parameters.
 The VCF files were then processed through biscuit vcf2bed, bedtools sort 
 (bedtools [37] version 2.29.2), and biscuit mergecg to create coordinate 
 sorted BED files containing CpG methylation beta value information




 biscuiteer
# Open biscuit files
 bisc <- readBiscuit(BEDfile = orig_bed, VCFfile = orig_vcf,
                    merged = FALSE)
bisc2 <- readBiscuit(BEDfile = shuf_bed, VCFfile = shuf_vcf,
                     merged = FALSE)

# Combine all 
comb <- unionize(bisc, bisc2)
