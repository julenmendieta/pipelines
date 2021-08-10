
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
                    make_option(c("-s", "--sample_file"), type="character", default=NULL, help="Path to edgeR sample file with same sampleName as the ones in input_file after removing string_start and string_end", metavar="path"),
                    make_option(c("-l", "--length_column"), type="integer", default=NULL, help="Index 1 based location of length column"),
                    make_option(c("-d", "--data_columns"), type="integer", default=NULL, help="Index 1 based location of the column from which we have count data"),
                    make_option(c("-S", "--string_start"), type="character", default=NULL, help="Left side string of data column names that you want to remove. In R columnames '/' changes to '.' so add last part of path to file with this changes"),
                    make_option(c("-E", "--string_end"), type="character", default=NULL, help="Right side string of data column names that you want to remove")
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
if (is.null(opt$sample_file)){
    print_help(opt_parser)
    stop("Sample file must be supplied.", call.=FALSE)
}
if (is.null(opt$length_column)){
    print_help(opt_parser)
    stop("Length column number must be supplied.", call.=FALSE)
}
if (is.null(opt$data_columns)){
    print_help(opt_parser)
    stop("Data starting column must be supplied", call.=FALSE)
}
if (is.null(opt$string_start)){
    startString <- "allProcessed.bamfiles.valid."
} else {
    startString <- opt$string_start
}
if (is.null(opt$string_end)){
    endString <- ".sort.rmdup.rmblackls.rmchr.bam"
} else {
    endString <- opt$string_end
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
sampleFile <- opt$sample_file

## Importing gene expression (features count) data
Dataset = read.table(fileIn, header = TRUE, sep = "\t",row.names = 1,
                          as.is = TRUE)

# remove extra info
lengths <- Dataset[,lenCol]
Dataset <- Dataset[,dataCol:ncol(Dataset)]

# reduce column names length
newColNames <- unlist(lapply(colnames(Dataset), function(x) 
                                strsplit(x, startString)[[1]][2]))
newColNames <- unlist(lapply(newColNames, function(x) 
                                strsplit(x, endString)[[1]][1]))
colnames(Dataset) <- newColNames

## Importing sample Info (drug treatment and cell lines info)
SampleInfo = read.table(sampleFile, header = TRUE, sep = "\t")
SampleInfo = as.matrix(SampleInfo)

# Using DGElist function from edgeR package so that it would be easy to do further analysis
DGE = edgeR::DGEList(counts = Dataset, lib.size = colSums(Dataset),
                     norm.factors = calcNormFactors(Dataset), samples = SampleInfo,group = NULL,
                     genes = NULL, remove.zeros = FALSE)
# more normalization factor value for less library size samples and viceversa was generated by
# calcNormFactors in the DGEList function

cpm = edgeR::cpm(DGE)

# get RPKM
#rpkm = edgeR::rpkm(DGE, gene.length=lengths)

write.table(data.frame("interval"=rownames(cpm),cpm), file=fileOut, quote=FALSE, sep="\t",
            row.names=FALSE)



