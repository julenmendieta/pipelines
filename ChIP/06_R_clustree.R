#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(reticulate)
    library(glue)
    library(clustree)
    library(optparse)

})


option_list <- list(make_option(c("-i", "--input_file"), type="character", 
                                default=NULL, 
                                help="Input CSV with cluster assignments in columns", 
                                metavar="path"),
                    make_option(c("-o", "--outfile"), type="character", default=NULL, 
                                help="Output image path", metavar="path"),
                    make_option(c("-c", "--nodeColor"), type="character", default=NULL, 
                                help="Column to use to color nodes"),
                    make_option(c("-n", "--nodeCount"), type="character", default=NULL, 
                                help="Column to use to show summed count")
                    )

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

opt$input_file
infile <- opt$input_file
OUT_PATH <- opt$outfile
nodeColor <- opt$nodeColor
nodeCount <- opt$nodeCount
if (is.null(nodeCount)) {
    node_label_aggr <- NULL
} else {
    node_label_aggr <- "sum"
}
setwd(dirname(OUT_PATH))
##################################################################
# Functions
##################################################################

set_fig_dimensions = function(num_clusterings){
    width = 10
    height = (0.6 * num_clusterings)
    
    if (height < 8){
        height = 8
    }
    
    png(width = width, height = height)
    options(repr.plot.width = width, repr.plot.height = height)
    
    return(list(width=width,height=height))
}

##################################################################
# Code
##################################################################
cluster_df <- read.csv(infile)

dims = set_fig_dimensions(num_clusterings = ncol(cluster_df))
# dims

# options(repr.plot.width = 10, repr.plot.height = 10)

if (is.null(nodeColor)) {
    g = clustree(
        x=cluster_df,
        prefix="K",
        node_label = nodeCount,
        node_label_aggr = node_label_aggr)
} else {
    g = clustree(
        x=cluster_df,
        prefix="K",
        node_colour=nodeColor,
        node_colour_aggr = "mean",
        node_label = nodeCount,
        node_label_aggr = node_label_aggr)
}
g


ggsave(
    OUT_PATH, 
    width = dims$width, 
    height = dims$height, 
    dpi = 100
)