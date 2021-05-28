#!/usr/bin/env Rscript
#===============================================================================
#' Author: Maria E. Calleja
#' Date: 2020/05
#' Recycled from MECC_sgRNA_analysis for run363-365
# David's comment: Julen, you only need to run the multimerge function to combine the individual tables (one table per sample) from the previous step into one table containing every sample from each experiment. After that you just need to split this large table into: a) table with target sgRNAs b) table with "noise" = sgRNAs that shouldn't be there. And calculate: 1) Read ounts/sample (Colsums removing the * entry that contains the unmapped reads) 2) boxplots (in log scale) and  3) density plots for each of them.
#===============================================================================
## PACKAGES.
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(tidyverse)
library(lazyeval)
library(dplyr)
library(tibble)
#===============================================================================
## GLOBAL VARIABLES.
## TO BE RUN @FANGORN with cluster connection
PARENT_DIR<-("/home/mecc/clustermecc")
PARENT_DIR<-("/home/bioinfo/Cluster")
PROJECT_DIR<-file.path(PARENT_DIR,"agoni_david_crispr")
RAW_DATA_DIR<-file.path(PROJECT_DIR,"data/edited_444")
SCRIPT_DIR<-file.path(PROJECT_DIR,"src")
TEMP_DIR<-file.path(PROJECT_DIR,"temp")
DL_DIR<-file.path(PROJECT_DIR,"data/edited_444")
Counts<-file.path(DL_DIR,"QC")
RSession<-file.path(PROJECT_DIR,"RSession_444")
#===============================================================================
######## FUNCTIONS. (the first is just an scaffold example from a known analysis)
#' @param x Character vector of length 1, f
#' @return Character vector of length 1, 
#' @export

multimerge <- function (mylist) {
  ## mimics a recursive merge or full outer join
  unames <- unique(unlist(lapply(mylist, rownames)))
  n <- length(unames)
  out <- lapply(mylist, function(df) {
    tmp <- matrix(nr = n, nc = ncol(df), dimnames = list(unames,colnames(df)))
    tmp[rownames(df), ] <- as.matrix(df)
    rm(df); gc()
    return(tmp)
  })
  stopifnot( all( sapply(out, function(x) identical(rownames(x), unames)) ) )
  bigout <- do.call(cbind, out)
  colnames(bigout) <- paste(rep(names(mylist), sapply(mylist, ncol)), unlist(sapply(mylist, colnames)), sep = "_")
  return(bigout)
}
#My though on this: is that tidyverse is here to help, but I need to re-cap. 

#===============================================================================
### MAIN.
################################################################################
## Step . data loading
################################################################################

#Set the directory to where your "extracted.fastq.bam.cnt.txt.final.txt" are located
setwd (DL_DIR)

#Read each of the files in a different table. See example below:
countreports <- dir(file.path(Counts), pattern = "*.indxstats")
#for (reporte in countreports){
#   t <- read.table(file.path(Counts, reporte), stringsAsFactors = FALSE,
#                   col.names= c("sgRNA", "length","counts", "unmapped"), row.names = 1) 
#                   #na.strings = c('', 'NA', '<NA>'),  blank.lines.skip = TRUE)
#   t <- t[,c("sgRNA", "counts", "unmapped")]
#   assign(reporte, t, envir = .GlobalEnv)
#   rm(t)
#}
countreports <- countreports [grep("_B_", countreports)]  
#dim(read.table(file.path(Counts,reporte)))

for (reporte in countreports){
  t <- read.table(file.path(Counts, reporte), stringsAsFactors = FALSE,
                  col.names= c("sgRNA","length","counts", "unmapped"), 
                  row.names = 1) 
  row.names(t)[nrow(t)] <- "unmapped"
  t[nrow(t),"counts"] <- t[length(rownames(t)), "unmapped"]
  t <- t[,"counts", drop = FALSE]
  assign(reporte, t, envir = .GlobalEnv)
  rm(t)
}

countable_list <- Filter(function(x) is(x, "data.frame"), mget(ls()))

#Use the multimerge function to merge them together into a single table.
table<- multimerge(countable_list)
colnames(table)
colnames(table) <- gsub (".indxstats_counts","", colnames(table))
# #Or else, #colnames(table) <- gsub ("_extracted.fastq.bam.cnt.txt.final.txt_COUNTS","", colnames(table))
# table <- table [-c(1),] #I think I ran it twice... for that reason I need to delete first row
# rownames(table)[1] <- "unmaped"
table <- data.frame(apply(table, 2, function(x) as.numeric(x)), row.names = rownames(table))
rm(list = ls(pattern = ".indxstats"))
table[is.na(table)] <- 0 # Replace NAs for zeros
write.csv(table, file.path(RSession,"CRISPR_RAW_B_TABLE.csv"), quote = FALSE, row.names = TRUE)

par(mar=c(14,5,1,1))
barplot(as.matrix(table[nrow(table),]), las=2, cex.names = 1.5, main = "Unaligned reads") 

pdf(file.path(RSession,"UnmappedReads_B.pdf"), width=11, height=8.5)
par(mar=c(14,5,1,1))
barplot(as.matrix(table[nrow(table),]), las=2, cex.names = 1.5, main = "Unaligned reads") 
dev.off()
par()
dev.off()

table<-table[-nrow(table),] # Remove row that contains unaligned reads
table<-table[apply(table,1,max)>1000,]
plot(density(unlist(log2(table+1))))

pdf(file.path(RSession,"Density_after1000_noumapped_B.pdf"), width=11, height=8.5)
plot(density(unlist(log2(table+1))))
dev.off()


## David's comment: Julen, finish here...the remaining is about extracting basic stats (you can code this yourself better) 1) Read ounts/sample (Colsums removing the * entry that contains the unmapped reads) 2) boxplots (in log scale) and  3) density plots for each of them.

##  David's comment:

################################################################################
## Step . data analysis
################################################################################

#Check the sgRNA count distribution:
pdf(file.path(RSession,"Summary_B.pdf"), width=11, height=8.5)
cols <- rainbow(10) # The number should be the same as the number of samples you have
par(mfrow = c(1, 1))
plot(NA, ylim = c(0, 1), xlim = c(2, 18), main="sgRNA Distribution before normalization") # You can change the xlim and ylim to fit your distribution
for (i in 1:ncol(table)){
  lines (density(log2(na.omit(table[,i]))+1), col=cols[i], lwd=2)
}
legend("topright", legend=colnames(table), lwd=2, col=cols)

par(mar=c(15, 4.1, 4.1, 2.1), cex=0.5)
boxplot(log2(table+1),las=2)
barplot(colSums(table), las=2, main="Total counts per sample")
dev.off()

write.csv(table, file.path(RSession,"CRISPR_RAW_B_FILTERED_TABLE.csv"), quote = FALSE, row.names = TRUE)
quantile(apply(table,1,max)) # Explore distribution

#Normalize data ##
##  David's comment: We are now using a better normalization method but this one is useful to have a quick look to the data. Anyway, this can be done by Ainhoa or me once you send us the raw tables.
factor.sizes<-1/(colSums(na.omit(table))/ mean(colSums(na.omit(table)))) # Calculates "factor.sizes"
table.norm<-t(t(table) * factor.sizes) # Applies normalisation
table.norm<-as.data.frame(table.norm)

pdf(file.path(RSession,"Summary_B_after_norm.pdf"), width=11, height=8.5)
cols <- rainbow(10)
par(mfrow = c(1, 1))
plot(NA, ylim = c(0, 1), xlim = c(10, 20), main="sgRNA Distribution after normalization")
for (i in 1:ncol(table.norm)){
  lines (density(log2(table.norm[,i])+1), col=cols[i], lwd=2)
}
legend("topright", legend=colnames(table.norm), lwd=2, col=cols)

par(mar=c(14, 4.1, 4.1, 2.1), cex=0.5)
boxplot(log2(table.norm+1),las=2)
barplot(colSums(table.norm), las=2, main="Total counts per sample")
dev.off()

write.csv(table.norm, file.path(RSession,"CRISPR_NORM_B_FILTERED_TABLE.csv"), quote = FALSE, row.names = TRUE)
