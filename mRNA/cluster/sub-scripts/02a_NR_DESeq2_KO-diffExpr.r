#!/usr/bin/env Rscript


################################################
################################################
## LOAD LIBRARIES                             ##
################################################
################################################

library(optparse)
library(DESeq2)
library(vsn)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(lattice)
library(BiocParallel)
library(limma)


################################################
################################################
## PARSE COMMAND-LINE PARAMETERS              ##
################################################
################################################

option_list <- list(make_option(c("-i", "--countTable"), type="character", 
                        default=NULL, 
                        help="Raw per gene read counts file.", metavar="path"),
                    make_option(c("-o", "--outdir"), type="character", 
                        default='./', help="Output directory", metavar="path"),
                    make_option(c("-p", "--outprefix"), type="character", 
                        default='differential', help="Output prefix", 
                        metavar="string"),
                    make_option(c("-c", "--cores"), type="integer", default=1, 
                        help="Number of cores", metavar="integer"),
                    make_option(c("-t", "--controls"), type="character", 
                        default='', help="comma separated control IDs", 
                        metavar="string"))

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

countTable <- opt$countTable
outdir <- opt$outdir
outprefix <- opt$outprefix
cores <- opt$cores
controls <- opt$controls

################################################
################################################
## READ IN COUNTS FILE                        ##
################################################
################################################

count.table <- read.delim(file=countTable,header=TRUE,check.names = FALSE)
rownames(count.table) <- count.table$gene_id
geneNames <- count.table['gene_name']
geneIds <- count.table['gene_id']
count.table = count.table[,colnames(count.table)[3:ncol(count.table)]]

# I'll have to fix this in the previous scripts in the future
# '-' characters have been changed to '.'
cnames <- colnames(count.table)
cnames <- gsub("\\.", "-", cnames)
colnames(count.table) <- cnames


################################################
################################################
## RUN DESEQ2                                 ##
################################################
################################################

if (file.exists(outdir) == FALSE) {
    dir.create(outdir,recursive=TRUE)
}
setwd(outdir)

samples.vec_ <- sort(colnames(count.table))
cells <- unique(sub("-.*", "\\1", samples.vec_))

for (cell in cells) {
    samples.vec <- grep(cell, samples.vec_, value=T)

    # JULEN: changed it to get as group the cell type and chip (or batch)
    #groups <- sub("_.*$", "", samples.vec)
    groups <- sub("^([^_]*_[^_]*).*", "\\1", samples.vec)
    batches <- gsub("^.*_", "", groups)
    groups <- sub("^([^_]*).*", "\\1", groups)
    labels <- gsub("_.*","",samples.vec)
    # Get the KO ID, we dont care about specific ko-guide IDs
    ko <- sapply(strsplit(labels,"-"), `[`, 2)
    


    # Ensure we only have on cell
    if (length(cell) != 1) {
        stop("Only same-cell comparisons are allowed.", call.=FALSE)
    }

    # convert controls string to list
    controls <- strsplit(controls, ',')[[1]]

    print(unique(groups))
    if (length(unique(groups)) == 1) {
        quit(save = "no", status = 0, runLast = FALSE)
    }


    counts <- count.table[,samples.vec,drop=FALSE]
    coldata <- data.frame(row.names=colnames(counts),
                condition=groups, batch=batches, labels=labels,
                ko=ko)
    coldata2 <- data.frame(row.names=colnames(counts),
                condition=groups, batch=batches, labels=labels,
                ko=ko)

    # Now we want to check if we have the same condition in more than one batch
    # and filter out one-batch comparisons and controls

    ## BATCH CORRECTION IS NOT IMPLEMENTED YET
    getBatchCorrected <- c()
    for (con in unique(coldata2[,"labels"])) {
        pos <- coldata2[,"labels"] == con
        if (length(unique(coldata2[pos, "batch"])) > 1) {
            getBatchCorrected <- c(getBatchCorrected, con)
            
        }
    }
    coldata2 <- coldata2[coldata2[,"labels"] %in% getBatchCorrected,]
    coldata2 = coldata2[!(coldata2[,'ko'] %in% controls),]
    coldata[coldata[,"ko"] %in% controls, "ko"] = 'Control'
    ko_list <- ko
    ko_list[ko %in% controls] = 'Control'

    # We iterate over labels to compare diff batches
    batchGroups <- c()
    i = 1
    for (la in unique(coldata2[,"labels"])) {
        pos <- coldata2[,'labels'] == la
        batchGroups[[i]] <- sort(unique(coldata2[pos,"batch"]))
        i <- i + 1
    }
    batchGroups <- unique(batchGroups)

    # If all same batch
    if (is.null(batchGroups)) {
        batchGroups <- unique(batches)
        coldata2 <- coldata[coldata['ko'] != 'Control',]
    } else {
        print("Batch correction not optimised")
        stopifnot(2 > 500)
    }

    for (bag in batchGroups) {
        bag <- sort(bag)
        compareID <- paste(c(cell, bag), collapse='-')
        DDSFile <- paste(outprefix,"_", compareID, ".dds.rld.RData",sep="")
        if (file.exists(DDSFile) == FALSE) {

            print(bag)
            coldataSel <- coldata2[coldata2[, "batch"] %in% bag, ]
            # Keep samples with data in all batches
            getBatchCorrected2 <- c()
            for (con in unique(coldataSel[,"labels"])) {
                pos <- coldata2[,"labels"] == con
                if (identical(sort(unique(coldata2[pos, "batch"])), bag)) {
                    getBatchCorrected2 <- c(getBatchCorrected2, con)
                    
                }
            }
            coldataSel <- coldataSel[coldataSel[, "labels"] %in% getBatchCorrected2, ]
            # Add controls
            coldataSel <- rbind(coldataSel, 
                coldata[(coldata[,"batch"] %in% bag) & (coldata[,"ko"] == "Control"), ])

            if (length(unique(coldataSel[, "batch"])) > 1) {
                dds <- DESeqDataSetFromMatrix(countData = round(counts[, rownames(coldataSel) ]), 
                            colData = coldataSel[, c("ko", "batch")], design = ~ batch + ko)
            } else {
                dds <- DESeqDataSetFromMatrix(countData = round(counts[, rownames(coldataSel) ]), 
                            colData = coldataSel[, c("ko", "batch")], design = ~ ko)
            }
            
            # Add gene Names
            mcols(dds) <- DataFrame(mcols(dds), geneNames)

            # Run diff analysis
            dds <- DESeq(dds, parallel=TRUE, BPPARAM=MulticoreParam(cores))
            # log2 transform count data minimizing differences between samples 
            #   for rows with small counts
            rld <- rlog(dds)
            vsd <- vst(dds, blind = FALSE)
            mat <- assay(vsd)
            if (length(unique(coldataSel[, "batch"])) > 1) {
                mat <- limma::removeBatchEffect(mat, vsd$batch)
                assay(vsd) <- mat
            } 
            save(dds,rld, vsd, file=DDSFile)
        } else {
            load(DDSFile)
            counts <- count.table[,samples.vec,drop=FALSE]
            coldataSel <- coldata2[coldata2[, "batch"] %in% bag, ]
            # Keep samples with data in both batches
            getBatchCorrected2 <- c()
            for (con in unique(coldataSel[,"labels"])) {
                pos <- coldata2[,"labels"] == con
                if (identical(sort(unique(coldata2[pos, "batch"])), bag)) {
                    getBatchCorrected2 <- c(getBatchCorrected2, con)
                    
                }
            }
            coldataSel <- coldataSel[coldataSel[, "labels"] %in% getBatchCorrected2, ]
            # Add controls
            coldataSel <- rbind(coldataSel, 
                coldata[(coldata[,"batch"] %in% bag) & (coldata[,"ko"] == "Control"), ])

        }


        ################################################
        ################################################
        ## PLOT QC                                    ##
        ################################################
        ################################################

        PlotFile <- paste(outprefix, "_", compareID, ".plots.pdf",sep="")
        if (file.exists(PlotFile) == FALSE) {
            pdf(file=PlotFile,onefile=TRUE,width=7,height=7)

            ## PCA
            pca.data <- DESeq2::plotPCA(rld,intgroup=c("ko"),returnData=TRUE)
            percentVar <- round(100 * attr(pca.data, "percentVar"))
            plot <- ggplot(pca.data, aes(PC1, PC2, color=ko)) +
                    geom_point(size=3) +
                    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
                    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
                    theme(panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(),
                        panel.border = element_rect(colour = "black", fill=NA, size=1))
            print(plot)

            ## WRITE PC1 vs PC2 VALUES TO FILE
            pca.vals <- pca.data[,1:2]
            colnames(pca.vals) <- paste(colnames(pca.vals),
                            paste(percentVar,'% variance',sep=""), sep=": ")
            pca.vals <- cbind(sample = rownames(pca.vals), pca.vals)
            write.table(pca.vals,file=paste(outprefix,"_", compareID,".pca.vals.txt",sep=""),
                        row.names=FALSE,col.names=TRUE,sep="\t",quote=TRUE)

            ## SAMPLE CORRELATION HEATMAP
            sampleDists <- dist(t(assay(rld)))
            sampleDistMatrix <- as.matrix(sampleDists)
            colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
            pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists,
                    clustering_distance_cols=sampleDists,col=colors)

            ## WRITE SAMPLE DISTANCES TO FILE
            write.table(cbind(sample = rownames(sampleDistMatrix), sampleDistMatrix),
                        file=paste(outprefix,"_", compareID,".sample.dists.txt",sep=""),
                        row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)

            dev.off()
        }

        if (length(unique(coldataSel[, "batch"])) > 1) {
            PlotFile <- paste(outprefix,"_", compareID, ".batchCorrected.plots.pdf",sep="")
            if (file.exists(PlotFile) == FALSE) {
                pdf(file=PlotFile,onefile=TRUE,width=7,height=7)

                ## PCA
                pca.data <- DESeq2::plotPCA(vsd,intgroup=c("ko"),returnData=TRUE)
                percentVar <- round(100 * attr(pca.data, "percentVar"))
                plot <- ggplot(pca.data, aes(PC1, PC2, color=ko)) +
                        geom_point(size=3) +
                        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
                        ylab(paste0("PC2: ",percentVar[2],"% variance")) +
                        theme(panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            panel.background = element_blank(),
                            panel.border = element_rect(colour = "black", fill=NA, size=1))
                print(plot)

                ## WRITE PC1 vs PC2 VALUES TO FILE
                pca.vals <- pca.data[,1:2]
                colnames(pca.vals) <- paste(colnames(pca.vals),paste(percentVar,'% variance',sep=""),
                                        sep=": ")
                pca.vals <- cbind(sample = rownames(pca.vals), pca.vals)
                write.table(pca.vals,file=paste(outprefix,"_", compareID, ".batchCorrected.pca.vals.txt",sep=""),
                        row.names=FALSE,col.names=TRUE,sep="\t",quote=TRUE)

                ## SAMPLE CORRELATION HEATMAP
                sampleDists <- dist(t(assay(vsd)))
                sampleDistMatrix <- as.matrix(sampleDists)
                colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
                pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists,
                        clustering_distance_cols=sampleDists,col=colors)

                ## WRITE SAMPLE DISTANCES TO FILE
                write.table(cbind(sample = rownames(sampleDistMatrix), sampleDistMatrix),
                            file=paste(outprefix,"_", compareID, ".batchCorrected.sample.dists.txt",sep=""),
                            row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)

                dev.off()
            }
        }
        ################################################
        ################################################
        ## SAVE SIZE FACTORS                          ##
        ################################################
        ################################################

        SizeFactorsDir <- "sizeFactors/"
        if (file.exists(SizeFactorsDir) == FALSE) {
            dir.create(SizeFactorsDir,recursive=TRUE)
        }

        NormFactorsFile <- paste(SizeFactorsDir,outprefix,"_", compareID, ".sizeFactors.RData",sep="")
        if (file.exists(NormFactorsFile) == FALSE) {
            normFactors <- sizeFactors(dds)
            save(normFactors,file=NormFactorsFile)

            for (name in names(sizeFactors(dds))) {
                sizeFactorFile <- paste(SizeFactorsDir,name,"_", compareID, ".sizeFactor.txt",sep="")
                if (file.exists(sizeFactorFile) == FALSE) {
                    write(as.numeric(sizeFactors(dds)[name]),file=sizeFactorFile)
                }
            }
        }

        ################################################
        ################################################
        ## WRITE LOG FILE                             ##
        ################################################
        ################################################

        LogFile <- paste(outprefix,"_", compareID, ".log",sep="")
        if (file.exists(LogFile) == FALSE) {
            cat("\nSamples =",samples.vec,"\n\n",file=LogFile,append=FALSE,sep=', ')
            cat("Groups =",groups,"\n\n",file=LogFile,append=TRUE,sep=', ')
            cat("Batches =",batches,"\n\n",file=LogFile,append=TRUE,sep=', ')
            cat("Dimensions of count matrix =",dim(counts),"\n\n",file=LogFile,append=TRUE,sep=' ')
            cat("\n",file=LogFile,append=TRUE,sep='')
        }


        ################################################
        ################################################
        ## LOOP THROUGH COMPARISONS                   ##
        ################################################
        ################################################

        raw.counts <- counts(dds,normalized=FALSE)
        colnames(raw.counts) <- paste(colnames(raw.counts),'raw',sep='.')
        pseudo.counts <- counts(dds,normalized=TRUE)
        colnames(pseudo.counts) <- paste(colnames(pseudo.counts),'pseudo',sep='.')
        if (length(unique(coldataSel[, "batch"])) > 1) {
            batchC.counts <- assay(vsd)
            colnames(batchC.counts) <- paste(colnames(batchC.counts),'unbatched',sep='.')
        }

        ResultsFile <- paste(outprefix, ".results_", compareID,".txt",sep="")

        addidx <- 1

        deseq2_results_list <- list()
        comparisons <- combn(unique(coldataSel[,'ko']),2)

        for (idx in 1:ncol(comparisons)) {

            # We only want to compare KO vs Controls
            if (sum(comparisons[,idx] %in% c("Control")) == 1) {

                control.group <- comparisons[1,idx]
                treat.group <- comparisons[2,idx]
                CompPrefix <- paste(paste0(cell, "-", control.group, "_", paste(bag, collapse='-')), 
                                    paste0(cell, "-", treat.group, "_", paste(bag, collapse='-')),sep="-vs-")
                cat("Saving results for ",CompPrefix," ...\n",sep="")

                CompOutDir <- paste(CompPrefix,'/',sep="")
                if (file.exists(CompOutDir) == FALSE) {
                    dir.create(CompOutDir,recursive=TRUE)
                }

                control.samples <- rownames(coldataSel[coldataSel['ko'] == control.group,])
                treat.samples <- rownames(coldataSel[coldataSel['ko'] == treat.group,]) 
                comp.samples <- c(control.samples,treat.samples)

                comp.results <- results(dds,contrast=c("ko",c(control.group,treat.group)))
                comp.df <- as.data.frame(comp.results)
                if (length(unique(coldataSel[, "batch"])) > 1) {
                    comp.table <- cbind(geneNames, as.data.frame(comp.df), 
                                raw.counts[,paste(comp.samples,'raw',sep='.')], 
                                pseudo.counts[,paste(comp.samples,'pseudo',sep='.')],
                                batchC.counts[,paste(comp.samples,'unbatched',sep='.')])
                } else {
                    comp.table <- cbind(geneNames, as.data.frame(comp.df), 
                                raw.counts[,paste(comp.samples,'raw',sep='.')], 
                                pseudo.counts[,paste(comp.samples,'pseudo',sep='.')])
                }


                ## WRITE RESULTS FILE
                CompResultsFile <- paste(CompOutDir,CompPrefix,".deseq2.results.txt",sep="")
                write.table(comp.table, file=CompResultsFile, col.names=TRUE, 
                            row.names=TRUE, sep='\t', quote=FALSE)

                ## FILTER RESULTS BY FDR & LOGFC AND WRITE RESULTS FILE
                if (length(comp.samples) > 2) {
                    pdf(file=paste(CompOutDir,CompPrefix,".deseq2.plots.pdf",sep=""),width=10,height=8)
                    for (MIN_FDR in c(0.01,0.05)) {

                        ## SUBSET RESULTS BY FDR
                        pass.fdr.table <- subset(comp.table, padj < MIN_FDR)
                        pass.fdr.up.table <- subset(comp.table, padj < MIN_FDR & log2FoldChange > 0)
                        pass.fdr.down.table <- subset(comp.table, padj < MIN_FDR & log2FoldChange < 0)

                        ## SUBSET RESULTS BY FDR AND LOGFC
                        pass.fdr.logFC.table <- subset(comp.table, padj < MIN_FDR & abs(log2FoldChange) >= 1)
                        pass.fdr.logFC.up.table <- subset(comp.table, padj < MIN_FDR & abs(log2FoldChange) >= 1 & log2FoldChange > 0)
                        pass.fdr.logFC.down.table <- subset(comp.table, padj < MIN_FDR & abs(log2FoldChange) >= 1 & log2FoldChange < 0)

                        ## WRITE RESULTS FILE
                        CompResultsFile <- paste(CompOutDir,CompPrefix,".deseq2.FDR",MIN_FDR,".results.txt",sep="")
                        CompBEDFile <- paste(CompOutDir,CompPrefix,".deseq2.FDR",MIN_FDR,".results.bed",sep="")
                        write.table(pass.fdr.table, file=CompResultsFile, col.names=TRUE, row.names=FALSE, sep='\t', quote=FALSE)

                        ## MA PLOT & VOLCANO PLOT
                        DESeq2::plotMA(comp.results, main=paste("MA plot FDR <= ",MIN_FDR,sep=""), ylim=c(-2,2),alpha=MIN_FDR)
                        plot(comp.table$log2FoldChange, -1*log10(comp.table$padj), col=ifelse(comp.table$padj<=MIN_FDR, "red", "black"), xlab="logFC", ylab="-1*log10(FDR)", main=paste("Volcano plot FDR <=",MIN_FDR,sep=" "), pch=20)

                        ## ADD COUNTS TO LOGFILE
                        cat(CompPrefix," genes with FDR <= ",MIN_FDR,": ",
                            nrow(pass.fdr.table)," (up=",
                            nrow(pass.fdr.up.table),", down=",
                            nrow(pass.fdr.down.table),")","\n",
                            file=LogFile,append=TRUE,sep="")
                        cat(CompPrefix," genes with FDR <= ",MIN_FDR," & FC > 2: ",
                            nrow(pass.fdr.logFC.table)," (up=",
                            nrow(pass.fdr.logFC.up.table),", down=",
                            nrow(pass.fdr.logFC.down.table),")","\n",
                            file=LogFile,append=TRUE,sep="")

                    }
                    cat("\n",file=LogFile,append=TRUE,sep="")
                    dev.off()
                }

                colnames(comp.df) <- paste(CompPrefix,".",colnames(comp.df),sep="")
                deseq2_results_list[[addidx]] <- comp.df
                addidx = addidx + 1

            }
            
            
        }
        ## WRITE RESULTS FROM ALL COMPARISONS TO FILE
        if (length(unique(coldataSel[, "batch"])) > 1) {
            deseq2_results_table <- cbind(geneNames,do.call(cbind, deseq2_results_list),
                                        raw.counts,pseudo.counts, batchC.counts)
        } else {
            deseq2_results_table <- cbind(geneNames,do.call(cbind, deseq2_results_list),
                                        raw.counts,pseudo.counts)
        }
        write.table(deseq2_results_table, file=ResultsFile, col.names=TRUE, 
                        row.names=TRUE, sep='\t', quote=FALSE)

    }
}

################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

RLogFile <- "R_sessionInfo.log"
if (file.exists(RLogFile) == FALSE) {
    sink(RLogFile)
    a <- sessionInfo()
    print(a)
    sink()
}

################################################
################################################
################################################
################################################
