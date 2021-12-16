#!/usr/bin/env Rscript

################################################
################################################
## LOAD LIBRARIES                             ##
################################################
################################################

library(optparse)
library(edgeR)

################################################
################################################
## PARSE COMMAND-LINE PARAMETERS              ##
################################################
################################################

option_list <- list(make_option(c("-i", "--input_file"), type="character", default=NULL, help="Path to featureCounts output ", metavar="path"),
                    make_option(c("-o", "--output_file"), type="character", default=NULL, help="Path to output CPM file", metavar="path"),
                    make_option(c("-l", "--length_column"), type="integer", default=NULL, help="Index 1 based location of length column"),
                    make_option(c("-d", "--data_columns"), type="integer", default=NULL, help="Index 1 based location of the column from which we have count data")
                    )

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input_file)){
    print_help(opt_parser)
    stop("Input file must be supplied.", call.=FALSE)
}
if (is.null(opt$output_file)){
    print_help(opt_parser)
    stop("Output file must be supplied.", call.=FALSE)
}
if (is.null(opt$length_column)){
    print_help(opt_parser)
    stop("Length column number must be supplied.", call.=FALSE)
}
if (is.null(opt$data_columns)){
    print_help(opt_parser)
    stop("Data starting column must be supplied", call.=FALSE)
}

################################################
################################################
## RUN CODE                                   ##
################################################
################################################

fileIn <- opt$input_file
fileOut <- opt$output_file
lenCol <- opt$length_column
dataCol <- opt$data_columns

## Importing gene expression (features count) data
Dataset = read.table(fileIn, header = TRUE, sep = "\t",row.names = 1,
                          as.is = TRUE, check.names = FALSE)

# remove extra info
lengths <- Dataset[,lenCol]
Dataset <- Dataset[,dataCol:ncol(Dataset)]

# reduce column names length following pipeline format
# (basename of first two '_' sections)
newColNames <- unlist(lapply(colnames(Dataset), function(x) 
                                tail(strsplit(x, '/')[[1]], n=1)))
newColNames <- unlist(lapply(newColNames, function(x)
                                strsplit(x, '\\.')[[1]][1]))
newColNames <- unlist(lapply(newColNames, function(x)
                                paste0(unlist(strsplit(x, '_')[[1]][1:2]), collapse='_')))
newColNames <- unlist(lapply(newColNames, function(x)
                                strsplit(x, '\\-')[[1]][1]))
                     
# if any, we remove NA introduced due to less than 3 sections
newColNames <- unlist(lapply(newColNames, function(x)
                                strsplit(x, '_NA')[[1]][1]))
colnames(Dataset) <- newColNames

## edgeR:: calcNormFactors
NormFactor <- calcNormFactors(object = Dataset, method = "TMM")

## if you prefer to use the DESeq2 strategy use method="RLE" instead

## raw library size:
LibSize <- colSums(Dataset)

## calculate size factors:
SizeFactors <- NormFactor * LibSize / 1000000

## Reciprocal, please read section below:   
SizeFactors.Reciprocal <- 1/SizeFactors
                             
# Store reciprocla values
write.table(SizeFactors.Reciprocal, file=fileOut, quote=FALSE, sep="\t",
            row.names=TRUE, col.names=FALSE)