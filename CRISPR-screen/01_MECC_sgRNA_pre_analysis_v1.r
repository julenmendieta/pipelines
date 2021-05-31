#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
#===============================================================================
#' Author: Maria E. Calleja
#' Date: 2020/05
#' Recycled from MECC_sgRNA_analysis for run363-365
# David's comment: Julen, you only need to run the multimerge function to combine 
# the individual tables (one table per sample) from the previous step into one 
# table containing every sample from each experiment. After that you just need 
# to split this large table into: a) table with target sgRNAs b) table with 
# "noise" = sgRNAs that shouldn't be there. And calculate: 1) Read counts/sample 
# (Colsums removing the * entry that contains the unmapped reads) 2) boxplots 
# (in log scale) and  3) density plots for each of them.
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
# How to run me
#Rscript --vanilla 01_MECC_sgRNA_pre_analysis_v1.r /Users/julen/Downloads/prueba/merge4-492 > r_out.txt
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>1) {
  stop("Too many arguments supplied (input file).n", call.=FALSE)
}
  
## GLOBAL VARIABLES.
## TO BE Modified
PROJECT_DIR <- args[1]
#PROJECT_DIR<-file.path("/Users/julen/Downloads/prueba/merge4-492")
#guidesFile <- "/Users/julen/Downloads/prueba/finalGuides.txt"
guidesFile <- "/home/jmendietaes/data/2021/CRISPR/finalGuides.txt"

# If the folder structure is ok this shouldn chance
Counts<-file.path(PROJECT_DIR,"idxstats")
RSession<-file.path(PROJECT_DIR,"RSession")
runName <- basename(PROJECT_DIR)

## Relevant info about the guides
allGuideCodes <- c("MGLibA", "R2.Br", "As", "B.Br", "TF1.Br", "Br", "R1.Br", "B", "TF2.Br")
guideSynonims <- vector(mode="list", length=length(allGuideCodes) + 5)
names(guideSynonims) <- c("TF1", "TF2", "R1", "R2", "BBr", allGuideCodes)
guideSynonims[[1]] <- "TF1.Br"; guideSynonims[[2]] <- "TF2.Br"
guideSynonims[[3]] <- "R1.Br"; guideSynonims[[4]] <- "R2.Br"
guideSynonims[[5]] <- "B.Br"

for (i in 1:length(allGuideCodes)) {
  guideSynonims[[i+5]] <- allGuideCodes[i]
}

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
setwd (PROJECT_DIR)

## Read each of the files in a different table. See example below:
# first we get all the files in folder and their associated libraries
countreports <- dir(file.path(Counts), pattern = "*.idxstats")
allLibs <- unique(sapply(strsplit(countreports, "_"), function(x) x[2]))
# sometimes we have files with more than one library separated by '-'
allLibs <- unique(unlist(strsplit(allLibs, "-")))
# all must be in guide or synonim list
if (! sum(allLibs %in% c(allGuideCodes, names(guideSynonims))) == length(allLibs)) {
  print('Some filetered guide libraries are not in the reference list')
  print(allLibs)
  # EXIT with error
  q("no", 1, FALSE)
}

for (glib_ in allLibs) {
  # get the real ID in case he file had a synonim
  glib = guideSynonims[[glib_]]
  
  print(paste0('---------------------- ', glib, ' ----------------------'))
  
  # Load the list that relates guides to each library
  guidesInfo <- read.table(guidesFile, stringsAsFactors = FALSE, header = TRUE) 
  guidesInfo <- guidesInfo[,c('ID', 'Library')]
  
  libPos <- grepl(glib, guidesInfo[,2])
  libGuides <- guidesInfo[libPos,1]
  
  # get the list of files for this library
  gcountreports <- countreports[grepl(glib_, countreports)]
  
  # read files in memory
  for (reporte in gcountreports){
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
  # keep only the data frames that contain the idxstats info
  countable_list <- countable_list[grepl("idxstats", names(countable_list))]
  
  
  # Use the multimerge function to merge them together into a single table.
  table<- multimerge(countable_list)
  colnames(table)
  colnames(table) <- gsub (".idxstats_counts","", colnames(table))
  # turn values to numbers
  table <- data.frame(apply(table, 2, function(x) as.numeric(x)), row.names = rownames(table))
  # delete lists form environment
  rm(list = ls(pattern = ".idxstats"))
  # Replace NAs for zeros
  table[is.na(table)] <- 0 
  # Write table with alignment info to all Libraries
  dir.create(RSession, showWarnings = FALSE)
  write.table(data.frame("ID"=rownames(table),table), paste0(RSession, "/", runName, "_", 
                                                             glib, "_crisprTable_raw_allLib.tsv"), 
              row.names=FALSE, quote=FALSE, sep='\t')
  
  
  # last line of table contains the unmaped guides
  pdf(paste0(RSession, "/", runName, "_", glib, "_stats.pdf"), width=11, height=8.5)
  
  
  #pdf(paste0(RSession, "/UnmappedReads_", glib, "_allLib.pdf"), width=11, height=8.5)
  par(mar=c(14,5,1,1))
  paste0("UnmappedReads_", glib, "_allLibrary")
  barplot(as.matrix(table[nrow(table),]), las=2, cex.names = 1.5, 
          main = paste0("UnmappedReads_", glib, "_allLibrary")) 
  #dev.off()
  #par(mar=c(14,5,1,1))
  #barplot(as.matrix(table[nrow(table),]), las=2, cex.names = 1.5, main = "Unaligned reads") 
  #par()
  #dev.off()
  # we now remove the unmaped reads list
  table<-table[-nrow(table),]
  
  
  ## Use valid guides info

  # First we split the data in valid and noise guides for this dataset
  # For that we get the library IDs that are contained in the selected library
  pos <- sapply(guidesInfo["Library"], function(x) grepl(glib, x))
  selectedIDs <- guidesInfo["ID"][pos]
  
  # Then we split the table as containing or not the valid guides 
  # and write the filtered elements
  validTable <- table[selectedIDs,]
  nonValidTable <- table[!(rownames(table) %in% selectedIDs),]
  
  
  write.table(data.frame("ID"=rownames(validTable),validTable), 
              paste0(RSession, "/", runName, "_", glib, 
                     "_crisprTable_raw_valid.tsv"), 
              row.names=FALSE, quote=FALSE, sep='\t')
  
  write.table(data.frame("ID"=rownames(nonValidTable),nonValidTable), 
              paste0(RSession, "/", runName, "_", glib, 
                     "_crisprTable_raw_nonValid.tsv"), 
              row.names=FALSE, quote=FALSE, sep='\t')

  
  
  # Remove guides with less than 1000 counts in all the samples
  #removed <- rownames(validTable[!apply(validTable,1,max)>1000,])
  #print(paste0(length(removed), ' guides removed for having less than 1000 counts in all samples'))
  #print(removed)
  #validTable<-validTable[apply(validTable,1,max)>1000,]
  #plot(density(unlist(log2(validTable+1))))
  
  
  #pdf(paste0(RSession, "/Density_", glib, "_valid.pdf"), width=11, height=8.5)
  plot(density(unlist(log2(validTable+1))), main = paste0("Density", glib, "_valid"))
  #dev.off()
  
  #pdf(paste0(RSession, "/Density_", glib, "_valid_after1000.pdf"), width=11, height=8.5)
  #plot(density(unlist(log2(validTable[apply(validTable,1,max)>1000,]+1))),
  #     main = paste0("Density", glib, "_valid_min1000"))
  #dev.off()
  
  #pdf(paste0(RSession, "/Density_", glib, "_NonValid.pdf"), width=11, height=8.5)
  plot(density(unlist(log2(nonValidTable+1))), main = paste0("Density", glib, "_nonValid"))
  #dev.off()
  
  #pdf(paste0(RSession, "/Density_", glib, "_NonValid_after1000.pdf"), width=11, height=8.5)
  #plot(density(unlist(log2(nonValidTable[apply(nonValidTable,1,max)>1000,]+1))),
  #     main = paste0("Density", glib, "_nonValid_min1000"))
  #dev.off()
  
  ## David's comment: Julen, finish here...the remaining is about extracting basic
  # stats (you can code this yourself better) 1) Read counts/sample (Colsums removing 
  # the * entry that contains the unmapped reads) 2) boxplots (in log scale) and  3) 
  # density plots for each of them.
  
  #Normalize data ##
  ##  David's comment: We are now using a better normalization method but this 
  # one is useful to have a quick look to the data. Anyway, this can be done by 
  # Ainhoa or me once you send us the raw tables.
  # Calculates "factor.sizes"
  factor.sizes<-1/(colSums(na.omit(validTable))/ mean(colSums(na.omit(validTable)))) 
  # Applies normalisation
  validTable.norm<-t(t(validTable) * factor.sizes) 
  validTable.norm<-as.data.frame(validTable.norm)
  
  write.table(data.frame("ID"=rownames(validTable.norm),validTable.norm), 
              paste0(RSession, "/", runName, "_", glib, 
                     "_crisprTable_normOld_valid.tsv"), 
              row.names=FALSE, quote=FALSE, sep='\t')
  
  ################################################################################
  ## Step . Basic stats
  ################################################################################
  ## For valid reads ##
  ## 1) Read counts/sample
  #Check the sgRNA count distribution:
  # RAW
  #pdf(paste0(RSession, "/readCount_", glib, "_valid_raw.pdf"), width=11, height=8.5)
  par(mar=c(15, 4.1, 4.1, 2.1), cex=0.5)
  barplot(colSums(validTable), las=2, main=paste0("Total counts per sample: ", glib))
  #dev.off()
  # NORM
  #pdf(paste0(RSession, "/CountBoxplot_", glib, "_valid_norm.pdf"), width=11, height=8.5)
  par(mar=c(15, 4.1, 4.1, 2.1), cex=0.5)
  barplot(colSums(validTable.norm), las=2, main=paste0("Total counts per sample: ", glib))
  #dev.off()
  
  
  ## 2) boxplots (in log scale)
  # RAW
  #pdf(paste0(RSession, "/CountBoxplot_", glib, "_valid_raw.pdf"), width=11, height=8.5)
  par(mar=c(15, 4.1, 4.1, 2.1), cex=0.5)
  boxplot(log2(validTable+1),las=2, main=paste0("Count boxplot per sample: ", glib))
  #dev.off()
  # NORM
  #pdf(paste0(RSession, "/CountBoxplot_", glib, "_valid_norm.pdf"), width=11, height=8.5)
  par(mar=c(15, 4.1, 4.1, 2.1), cex=0.5)
  boxplot(log2(validTable.norm+1),las=2, main=paste0("Normalised count boxplot per sample: ", glib))
  #dev.off()
  
  
  ## 3) density plots for each of them
  # RAW
  #pdf(paste0(RSession, "/CountDistrib_", glib, "_valid_raw.pdf"), width=11, height=8.5)
  cols <- rainbow(dim(validTable)[2]) 
  par(mfrow = c(1, 1))
  # You can change the xlim and ylim to fit your distribution
  plot(NA, ylim = c(0, 1), xlim = c(2, 18), 
       main=paste0("sgRNA Distribution before normalization: ", glib)) 
  for (i in 1:ncol(validTable)){
    lines (density(log2(na.omit(validTable[,i]))+1), col=cols[i], lwd=2)
  }
  legend("topright", legend=colnames(validTable), lwd=2, col=cols)
  #dev.off()
  # NORM
  #pdf(paste0(RSession, "/CountDistrib_", glib, "_valid_norm.pdf"), width=11, height=8.5)
  cols <- rainbow(dim(validTable.norm)[2]) 
  par(mfrow = c(1, 1))
  # You can change the xlim and ylim to fit your distribution
  plot(NA, ylim = c(0, 1), xlim = c(2, 18), 
       main=paste0("sgRNA Distribution after normalization: ", glib)) 
  for (i in 1:ncol(validTable.norm)){
    lines (density(log2(na.omit(validTable.norm[,i]))+1), col=cols[i], lwd=2)
  }
  legend("topright", legend=colnames(validTable.norm), lwd=2, col=cols)
  dev.off()
  
  
  
  # read percentyles valid vs non valid
  print('Max read count per guide percentyle in valid')
  print(quantile(apply(validTable,1,max))) # Explore distribution
  print('Max read count per guide percentyle in non-valid')
  print(quantile(apply(nonValidTable,1,max))) # Explore distribution
  
  # show the non valids with more than 10000 counts
  #p95 <- quantile(apply(nonValidTable,1,max), probs=c(0.95))
  print('The following non-valid guides have more thn 10k counts')
  print(nonValidTable[rowSums(nonValidTable > 10000) > 0, ])
  
}
