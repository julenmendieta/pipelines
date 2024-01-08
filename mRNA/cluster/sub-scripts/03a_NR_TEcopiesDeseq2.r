#!/usr/bin/env Rscript

################################################
################################################
## LOAD LIBRARIES                             ##
################################################
################################################

suppressPackageStartupMessages(library(DESeq2))
library(BiocParallel)
library(optparse)


################################################
################################################
## PARSE COMMAND-LINE PARAMETERS              ##
################################################
################################################

option_list <- list(make_option(c("-i", "--input_file"), type="character", default=NULL, help="htseq-count file", metavar="path"),
                    make_option(c("-o", "--outfile"), type="character", default='./TEcopies_Deseq2.txt', help="Output file with path", metavar="path"),
                    make_option(c("-c", "--cores"), type="integer", default=1, help="Number of cores", metavar="integer"),
                    make_option(c("-t", "--controls"), type="character", default='NTC', help="comma separated control IDs", metavar="string"))

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)


infile <- opt$input_file
outfile <- opt$outfile
controls <- opt$controls
cores <- opt$cores
################################################
################################################
## Run                                        ##
################################################
################################################

# Load data table
data <- read.table(infile,header=T,row.names=1,
                        check.names = FALSE)

min_read <- 1
data <- data[apply(data,1,function(x){max(x)}) > min_read,]

# convert controls string to list
controls <- strsplit(controls, ',')[[1]]
# Get sample groups
samples.vec <- sort(colnames(data))
# JULEN: changed it to get as group the cell type and chip (or batch)
#groups <- sub("_.*$", "", samples.vec)
groups <- sub("^([^_]*_[^_]*[^_]*).*", "\\1", samples.vec)
batches <- gsub("^.*_", "", groups)
groups <- sub("^([^-]*-[^_]*).*", "\\1", groups)
cells <- gsub("-.*", "", groups)


for (cell in unique(cells)) {
    cell.vec <- grep(paste0(cell, '-'), samples.vec, value=T)
    controlsPos <- unlist(lapply(controls, function(x) {grep(x, cell.vec, value=F)}))
    controlsCell <- cell.vec[controlsPos]
    kosCell <- cell.vec[-controlsPos]
    kos <- unlist(lapply(kosCell, function(x) {unlist(strsplit(x, '-'))[2]}))
    ukos <- unique(kos)
    
    for (ko in ukos) {
        id1 <- paste0(cell, '_', ko, '-vs-Control')
        TGroup <- grep(ko, kosCell, value=T)
        CGroup <- controlsCell
        groups_ <- c(TGroup, CGroup)
        data2 <- data[groups_]
        sampleInfo <- data.frame(groups_,row.names=colnames(data2))
        sampleInfo['ko'] <- factor(c(rep(ko,2),rep("Control",2)))
        
        # Get FC
        dds <- DESeqDataSetFromMatrix(countData = data2, colData = sampleInfo, design = ~ ko)
        dds$ko = relevel(dds$ko,ref="Control")
        dds <- DESeq(dds, parallel=TRUE, BPPARAM=MulticoreParam(cores))
        res <- results(dds)

        # rename and store new columns
        res <- as.data.frame(res)
        colnames(res) <- paste(id1, colnames(res), sep='.')
        data <- merge(data, as.data.frame(res), 'row.names')
        rownames(data) <- data[,'Row.names']
        data$Row.names <- NULL

    }
}


write.table(data, file=outfile, sep="\t",quote=F)