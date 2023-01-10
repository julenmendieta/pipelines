# Load libraries
library(csaw)
# http://bioconductor.org/books/3.15/csawBook/
library(edgeR)
library(BiocParallel)

########## INPUT VARIABLES #########
## Input variables
bamsPath <- "/scratch/julen/ChIP/bamFiles/main"
controlBamsPath <- "/scratch/julen/ChIP/bamFiles/IgG"
outPath <- "/scratch/julen/ChIP/allData/04_subsamplingNoIgG/outdata/csaw"

allowedCompare = list(list('DM', 'Mye'),
                      list('DM', 'Mye', 'Mono', 'MEP', 'GMP', 'Bcell', 'Ery')
                     )


cellOrder <- c('DM', 'Mye', 'LSKvitro', 'GMP', 'GMPvitro', 'Mono', 'MEP', 'MEPvitro', 'Ery', 
             'Bcell')
# Number of CPU's we will allow to use for this task
nCPU <- 24

# EL WINDOW SIZE VARIA CON EL TIPO DE PICO, IGUAL TENGO QUE MIRAR LA LONGITUD DE PICOS ESTIMADA DE MACS
# Size of the genomic bins (default 10 bp)
# Spacing between bins (middle is not check, but Even if spacing is 
#50, the reads are longer than that, so reads lying in the interval 
#between adjacent windows will still be counted into several windows)
# For histone marks, widths of at least 150 bp are recommended (Humburg et al. 2011). 
#This corresponds to the length of DNA wrapped up in each nucleosome
# The choice of window size controls the compromise between spatial resolution and count size. 
#Larger windows will yield larger read counts that can provide more power for DB detection.
# We suggest performing an initial DB analysis with small windows to maintain spatial resolution. 
#The widths of the final merged regions (see Section~) can provide an indication of the appropriate 
#window size. Alternatively, the analysis can be repeated with a series of larger windows, and the 
#results combined (see Section~). This examines a spread of resolutions for more comprehensive 
#detection of DB regions.
windowSize <- 10
# the start positions for adjacent windows are spacing bp apart.
spacing <- 50

# size of the bins we check to optimally filter out genomic locations
filterBinSize <- c(2000, 5000, 8000, 11000)

# Filtering is best done with globalEnrichment or localEnrichment (similar to MACS)
filterWay = 'globalEnrichment'
# This will assume that the controls are IgG, and that the experiments are done by cell type 
filterIgG <- TRUE
# AÑADIR OPCION D EILTRAR QUEDANDONOS CON LOS BINS QUE CAEN EN LOS CONSENSUS PEAKS
# have to convert it to GRanges object
## GRanges object with 6 ranges and 1 metadata column:
##             seqnames              ranges strand |     gene_id
##                <Rle>           <IRanges>  <Rle> | <character>
##   100009600     chr9   21062393-21076096      - |   100009600
##   100009609     chr7   84935565-84967115      - |   100009609
##   100009614    chr10   77708457-77712009      + |   100009614
# Path to the onsensus coordiantes file that should be transformed to Granges object
consensusCoord <- ''

# FC limit to mantain a bin when compared to IgG
IgG.fc <- 2

# minimum p-value for obtaining FC values
minFCpval <- 0.01

##############################################

## get all files and chip experiments in there
allFiles <- list.files(bamsPath)
allFiles <- allFiles[grep("*bam$", allFiles)]

if (filterIgG == TRUE) {
    controlFiles <- list.files(controlBamsPath)
    controlFiles <- grep("*bam$", controlFiles, value = TRUE)
    
    extraControlFiles <- list.files(paste(controlBamsPath, "cellMergeIgG", sep='/'))
    extraControlFiles <- grep("*bam$", extraControlFiles, value = TRUE)
}

# Look together only groups of cells we stated that can be compared
for (compare in allowedCompare) {
    print(compare)
    # the plus symbol in names is a problem, so we add some special characters on top
    comparePlusSafe <- gsub("\\+", "\\\\+", unlist(compare))
    chipFilesComp <- unlist(lapply(comparePlusSafe, function(x) grep(paste0(x, "_", collapse='_'), allFiles, value=TRUE)))
    cellGroup <- paste(unlist(compare), collapse='-')
                                                                  
    ## define CPU's to use
    multicoreParam <- MulticoreParam(workers = min(length(chipFilesComp), 64))
                          
    # make info dataframe
    Name <- unlist(lapply(chipFilesComp, function(x)
                            strsplit(x, '\\.')[[1]][1]))
    Name <- unlist(lapply(Name, function(x) {
                            x <- unlist(strsplit(x, '_'))[1:3]
                            x <- x[!is.na(x)]
                            paste(x, collapse='_')
                            }
        )
      )
    Description <- unlist(lapply(chipFilesComp, function(x)
                            strsplit(x, '_')[[1]][1]))
    Path <- paste(bamsPath, chipFilesComp, sep='/')
    infoData <- data.frame(Name, Description, Path)
    # sort table by cel order
    infoData <- infoData[order(match(infoData[, "Description"],  cellOrder)), ]               
    # With this variable we track the number of non IgG files
    chipPos <- nrow(infoData)
                                 
    ## get IgG control path
    if (filterIgG == TRUE) {
        for (cP in comparePlusSafe) {
            control <- unlist(lapply(cP, function(x) grep(x, controlFiles, value=TRUE)))
            # if no control files for this cell
            if (length(control) == 0) {
                if (! ("cellMergeIgG" %in% infoData[, "Description"])) {
                    control <- paste(controlBamsPath, "cellMergeIgG", extraControlFiles, sep='/')
                }
                
            } else {
                control <- paste(controlBamsPath, control, sep='/')
            }
            infoData <- rbind(infoData, c('IgG', cP, control))
            
        }
        
    }
    
    # Open blacklist regions (already filtered in BAM)
    #df1 <- read.table(blackList, sep="\t")
    # split and convert per region
    #blacklist_discard <- 
    #  lapply(split(df1, df1$V4), function(i){
    #    GRanges(seqnames = i$V1,
    #            ranges = IRanges(start = i$V2,
    #                             end = i$V3,
    #                             names = i$V4))})
                                     
    # We set the minimum mapping quality score to 
    # 10 to remove poorly or non-uniquely aligned reads
    # pe indicates that we take both pair-end reads
    param <- readParam(minq=10, pe='both')

    # This step fails with many reads, and is not needed for paired-end
    ## Computing the average fragment length of only ChIP experiments
    # For that we obtain the strand cross-correlaton 
    #x <- correlateReads(infoData[1:chipPos, 'Path'], param=reform(param, dedup=TRUE),
    #                   BPPARAM=multicoreParam)
    #frag.len <- maximizeCcf(x)
    #plot(1:length(x)-1, x, xlab="Delay (bp)", ylab="CCF", type="l")
    #abline(v=frag.len, col="red")
    #text(x=frag.len, y=min(x), paste(frag.len, "bp"), pos=4, col="red")

    ## Count reads by windows
    win.data <- windowCounts(infoData$Path, param=param, 
                                width=windowSize,  # bin size
                                #ext=frag.len,      # average length of sequence fragments (Only single-end)
                                spacing=spacing,   # distance between consecutive windows
                                filter=10,         # min bin reads (from all bamfile) to exclude
                                BPPARAM=multicoreParam)  # CPUs to use       

    #* win.data = RangedSummarizedExperiment object where the matrix of counts is stored as 
    #the first assay. Each row corresponds to a genomic window while each column corresponds 
    #to a library. The coordinates of each window are stored in the rowRanges. The total number 
    #of reads in each library (also referred to as the library size) is stored as totals in the colData

    # If you would want to get binned data from whole genome (no spacing)
    # demo <- windowCounts(infoData$Path, width=1000, bin=TRUE, param=param)
    # just think about the ram, hence enough width size. You can alternatively set 
    #filter=0 and allow spacing

    ## Filtering low abundance regions (we do it by global enrichment)
    # Also posible by count size, proportion, and enrichment
    if (filterWay == 'globalEnrichment') {
        # find the best filterBinSize
        print("Looking for best window width for filtering")
        TMMwindow <- filterBinSize[1]
        demo <- windowCounts(infoData[1:chipPos, 'Path'], bin=TRUE, width=TMMwindow, 
                             param=param, BPPARAM=multicoreParam)
        medDev <- median(abs(1- normFactors(demo, se.out=FALSE)))
        print(TMMwindow);print(medDev)

        for (i in filterBinSize[2:length(filterBinSize)]) {
            demo_ <- windowCounts(infoData[1:chipPos, 'Path'], bin=TRUE, width=i, 
                                  param=param, BPPARAM=multicoreParam)
            medDev_ <- median(abs(1- normFactors(demo_, se.out=FALSE)))

            print(i);print(medDev_)
            if (medDev_ < medDev) {
                demo <- demo_
                medDev <- medDev_
                TMMwindow <- i
            }
        }

        bins <- demo
        filter.stat <- filterWindowsGlobal(win.data[,1:chipPos], bins)

        min.fc <- 3
        hist(filter.stat$filter, main="", breaks=50,
                    xlab="Background abundance (log2-CPM)")
                    abline(v=log2(min.fc), col="red")

        
        keep <- filter.stat$filter > log2(min.fc)
        print(paste0("Summary of filtered ", TMMwindow, "bins"))

    # This is more similar to MACS
    } else if (filterWay == 'localEnrichment') {
        surrounds <- 2000
        neighbor <- suppressWarnings(resize(rowRanges(win.data[,1:chipPos]), surrounds, fix="center"))
        wider <- regionCounts(infoData[1:chipPos, 'Path'], regions=neighbor, 
                              #ext=frag.len, # Onlye single-end
                              param=param)
        filter.stat <- filterWindowsLocal(win.data[,1:chipPos], wider)

        hist(filter.stat$filter, main="", breaks=50,
                    xlab="Background abundance (log2-CPM)")
                    abline(v=log2(min.fc), col="red")

        min.fc <- 3
        keep <- filter.stat$filter > log2(min.fc)
        print(paste0("Summary of filtered ", TMMwindow, "bins"))

    } else if (filterWay == 'consensusPeaks') {
        print("No Statistic bin filtering done")
        print("Filtered by consensus peaks instead")
        suppressWarnings(keep <- overlapsAny(rowRanges(win.data[,1:chipPos]), consensusCoord))

    } else {
        print('Wrong filtering way')
        ERROR
    }

    summary(keep)
                                 
    # If we also want to filter by IgG signal
    if (filterIgG == TRUE) {
        print(paste0("Filtering by IgG FC > ", IgG.fc))
        
        # If we didn't run TMM filtering we need to set a window size
        if (exists("TMMwindow") == FALSE) {
            TMMwindow <- filterBinSize[1]
        }
            
            
        chip <- win.data[,1:chipPos] # All ChIP libraries
        control <- win.data[,chipPos+1] # All control libraries

        in.binned <- windowCounts(infoData$Path, bin=TRUE, width=TMMwindow, 
                                     param=param, BPPARAM=multicoreParam)

        chip.binned <- in.binned[,1:chipPos]
        control.binned <- in.binned[,chipPos+1]
        scale.info <- scaleControlFilter(chip.binned, control.binned)

        # filter with a prior count of 5 so small read values do not make huge FC
        filterControl.stat <- filterWindowsControl(chip, control, 
                                                    prior.count=5, scale.info=scale.info)
        # keep these regions and the ones that passed the enrichment filter
        keep <- keep & (filterControl.stat$filter > log2(IgG.fc))

    }
                                                   

    summary(keep)
    filtered.data <- win.data[keep,1:chipPos]
    #* How to acces data inside filtered.data
    #* matrix of counts
    #* head(assay(filtered.data))
    #* genomic coordinates
    #* head(rowRanges(filtered.data))

    ## NORMALISE BY COMPOSITION BIASES
    adjc <- calculateCPM(filtered.data, use.norm.factors=TRUE)

                                     
    ## write data file to tsv
    # First get data in dataframe format
    read.counts <- as.data.frame(assay(filtered.data))
    colnames(read.counts) <- paste0(infoData[1:chipPos, 'Name'], '.rcount')
    cpm.counts <- as.data.frame(adjc)
    colnames(cpm.counts) <- paste0(infoData[1:chipPos, 'Name'], '.cpm')
    read.counts <- cbind(as.data.frame(rowRanges(filtered.data)), read.counts, cpm.counts)
    read.counts$strand <- NULL
    # rename columns for coherence with consensus peaks and add interval_id column
    colnames(read.counts)[which(names(read.counts) == "seqnames")] <- "chr"
    read.counts["interval_id"] <- paste0("Interval_" , 1:nrow(read.counts))
    columns <- colnames(read.counts)[1:ncol(read.counts)-1]
    columns <- c(columns[1:3], 'interval_id', columns[4:length(columns)])
    read.counts <- read.counts[,columns]
                                     

    ## Statistical modelling
    # get normalisation offsets
    filtered.data <- normFactors(bins, se.out=filtered.data)
    normfacs <- normFactors(bins, se.out=FALSE)

    #y_base <- asDGEList(filtered.data, norm.factors=normfacs)
    #summary(y_base)


    # Construct design matrix for experimental design
    # TO BE CHANGED. in MY CASE WILL GO BY CHIP, IN THIS WAY WILL ALSO WORK
    #WHEN I ADD THE REPLICATES
    # https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/designmatrices.html#design-matrices-with-and-without-intercept-term
    celltype <- infoData$Description[1:chipPos]

    #celltype <- factor(celltype)
    #design <- model.matrix(~0+celltype)
    #colnames(design) <- levels(celltype)

    # we want to only compare same ChIP situations
    chips <- unlist(lapply(infoData$Name[1:chipPos], function(x)
                                strsplit(x, '_')[[1]][2]))
    chips <- unlist(lapply(chips, function(x)
                                strsplit(x, '-')[[1]][1])) 
    #chips <- factor(chips)
    #for (coln in colnames(design)) {
    #    design[,coln] <- design[,coln] * as.numeric(chips)
    #}

    #design

    ############################                             
    ## DIFFERENTIAL ANALYSIS
    ############################
    print("Differential analysis")
    # we get FC values when we have a ChIP more than once
    more1Chip <- as.character(as.data.frame(table(chips))[as.data.frame(table(chips))$Freq > 1, "chips"])
    for (chipComp in more1Chip) {
        print(chipComp)
        pos <- grep(chipComp, infoData[,'Name'])
        # if more than one cell
        if (length(unique(infoData[pos,'Description'])) > 1) {
            celltypeAll <- unique(infoData$Description[pos])

            cellPairs <- combn(celltypeAll, 2, simplify = F)
            for (pair in cellPairs) {
                print(unlist(pair))

                posCell <- grep(paste0(paste(unlist(pair), collapse='_|'), '_'), infoData[,'Name'])
                posCell <- posCell[posCell %in% pos]

                celltype <- infoData$Description[posCell]


                y_base <- asDGEList(filtered.data[,posCell], norm.factors=normfacs[posCell])

                # get design table for this comparison
                # in design the element with a 1 in celltype is divided by the element
                #with a 0. Due to that we reverse the order of cells in factors, so that
                #first cell is always divided by second
                celltype <- factor(rev(celltype))
                design <- model.matrix(~celltype)
                colnames(design) <- c('Intercept', 'celltype')

                # get differential data
                # supplying a “reasonable” value for the NB dispersion during GLM fitting (e.g., 0.05 - 0.1, 
                #based on past experience). DB windows are then identified using the likelihood ratio test.
                fit.norep <- glmFit(y_base, design, dispersion=0.075)
                results.norep <- glmLRT(fit.norep, contrast=c(0, 1))
                # add FC to common table
                compareId <- paste0(paste(paste(infoData[posCell,'Description'], chipComp, sep='_'), 
                                      collapse='-VS-'), '.log2FC')
                read.counts[compareId] <- results.norep$table$logFC
                # remove comparisons that were bellow p-value thresshold
                read.counts[results.norep$table$PValue >= minFCpval, compareId] <- NA

            }
        }
    }

                                 
    if (filterIgG == TRUE) {
        outfile <- paste0(outPath, '/binnedPeaks/allChIPCounts_', cellGroup, '_IgG-flt.tsv')
    } else {
        outfile <- paste0(outPath, '/binnedPeaks/allChIPCounts_', cellGroup, '.tsv')
    }
    write.table(read.counts, file=outfile, 
                    row.names=FALSE, quote=FALSE, sep="\t")
                                 
                                     
}
