library(methylKit)

# Input files can be gziped
inpath <- "/scratch/julen/DNAme/allData/01_firstTry/otdata/methylDackel"
species = 'mm10'
controlCell <- 'Mye'



infiles <- list.files(inpath)
for (methType in c('CpG')) {
    infiles_ <- grep('CpG', infiles, value=T)
    Ids <- (unlist(lapply(infiles_, function(x) 
            {paste(strsplit(x, '_' )[[1]][1:3], 
                    collapse='_')})))
    files <- file.path(inpath, infiles_)

    # get cell types 
    cells <- unlist(lapply(infiles_, function(x) 
            {strsplit(x, '_' )[[1]][1]}))
    batches <- unlist(lapply(infiles_, function(x) 
            {strsplit(x, '_' )[[1]][2]}))
    # if we have a preferred control, reorder all to have it first
    if (controlCell %in% cells) {
        Ids <- c(Ids[cells == controlCell], Ids[cells != controlCell]) 
        files <- c(files[cells == controlCell], files[cells != controlCell]) 
        batches <- c(batches[cells == controlCell], batches[cells != controlCell]) 
        cells <- c(cells[cells == controlCell], cells[cells != controlCell]) 
    }
    # and convert to numeric
    factorCells <- factor(cells, levels=unique(cells))
    factorCells <- as.numeric(factorCells) - 1


    # read the files to a methylRawList object: myobj
    myobj=methRead(as.list(files)  ,
            sample.id=as.list(Ids),
            assembly=species,
            treatment=factorCells,
            context=methType
            )

    ## Check data distribution
    # Normally, we should see bimodal distributions. Strong deviations from the 
    # bimodality may be due to poor experimental quality, such as problems with 
    # bisulfite treatment.
    for (nc in 1:length(cells)) {
        getMethylationStats(myobj[[nc]],plot=TRUE,both.strands=FALSE)
    }

    ## Get coverage
    # The bases with unusually high coverage are usually alarming. It might 
    # indicate a PCR bias issue in the experimental procedure
    for (nc in 1:length(cells)) {
        getCoverageStats(myobj[[nc]],plot=TRUE,both.strands=FALSE)
    }

    ## Filter samples based on coverage
    # They usually filter by minimum coverage of 10, that is reasonable, and 
    # then by percentyle 99.9. I checked the regions filtered by percentyle,
    # and some were increasing transcription after Trim28 KO, so wont use apply
    # this filter
    filtered.myobj=filterByCoverage(myobj,lo.count=10,lo.perc=NULL,
                                      hi.count=NULL,hi.perc=NULL)

    ## Create unified table with CpGs covered in at least 1 sample per group
    # Setting destrand=TRUE (the default is FALSE) will merge reads on both 
    # strands of a CpG dinucleotide (only advised when looking at CpG 
    # methylation (for CpH methylation this will cause wrong results in 
    # subsequent analyses))
    # min.per.group = an integer denoting minimum number of samples per 
    # replicate needed to cover a region/base. missing data for uncovered 
    # bases/regions will appear as NAs
    if (methType == "CpG") {
        destrand=FALSE
    } else {
        destrand=TRUE
    }
    meth=methylKit::unite(filtered.myobj, 
                            destrand=destrand, min.per.group=2L)

    ## Filter CpGs
    # pm=percMethylation(meth) # get percent methylation matrix
    # mds=matrixStats::rowSds(pm) # calculate standard deviation of CpGs
    # head(meth[mds>20,])
    # hist(mds,col="cornflowerblue",xlab="Std. dev. per CpG")

    # CC
    getCorrelation(meth,plot=TRUE)

    ## Clustering
    par(mar=c(8,4,4,1)+.1)
    clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)

    # PCA
    PCASamples(meth, screeplot=TRUE)
    PCASamples(meth, adj.lim=c(1,1))

    ## BATCH correction (TO BE TESTED)
    # All of these functions described in this section work on a methylBase 
    # object that does not have missing values (that means all bases in 
    # methylBase object SHOULD HAVE COVERAGE IN ALL SAMPLES)
    # Look for posible batch effects
    # sampleAnnotation=data.frame(batch_id=batches
    #                         #,age=c(19,34,23,40)
    #                         )
    # as=assocComp(mBase=meth,sampleAnnotation)
    # # construct a new object by removing the first pricipal component
    # # from percent methylation value matrix
    # newObj=removeComp(meth,comp=1)

    ## Tiling windows (TO BE TESTED)
    # For some situations, it might be desirable to summarize methylation 
    # information over tiling windows rather than doing base-pair resolution 
    # analysis.
    # The function below tiles the genome with windows of 1000bp length and 
    # 1000bp step-size and summarizes the methylation information on those 
    # tiles. In this case, it returns a methylRawList object which can be fed 
    # into unite and calculateDiffMeth functions consecutively to get 
    # differentially methylated regions
    # In this scenario one might want to set the initial per base coverage 
    # threshold of 10 to a lower value and then filter based on the number 
    # of bases (cytosines) per region
    if (2 == 3) {
        myobj_lowCov = methRead(as.list(files) ,
            sample.id=as.list(Ids),
            assembly=species,
            treatment=factorCells,
            context=methType,
            mincov = 3
            )
        tiles = tileMethylCounts(myobj_lowCov, win.size=1000, step.size=1000, 
                                cov.bases = 10)
        meth_lowCov=methylKit::unite(tiles, 
                                destrand=destrand, min.per.group=2L)
        myDiff_lowCov=calculateDiffMeth(meth_lowCov)
        diffMethPerChr(myDiff_lowCov,plot=FALSE,qvalue.cutoff=0.05, meth.cutoff=15,
                        type="all")
        getMethylDiff(myDiff_lowCov,difference=15,qvalue=0.05)
    }
    
    
    ## Differentila analysis
    myDiff=calculateDiffMeth(meth)
    # get all differentially methylated bases
    # If you specify type="hyper" or type="hypo" options, you will get 
    # hyper-methylated or hypo-methylated regions/bases
    myDiff25p=getMethylDiff(myDiff,difference=15,qvalue=0.05)
    diffMethPerChr(myDiff,plot=FALSE,qvalue.cutoff=0.05, meth.cutoff=15,
                    type="all")

    # Logistic regression based test
    my.diffMeth<-calculateDiffMeth(meth,
                                overdispersion="MN",test="Chisq",mc.cores=1)
    myDiff25p_2=getMethylDiff(my.diffMeth,difference=15,qvalue=0.05)
    diffMethPerChr(my.diffMeth,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=15,
                    type="all")

    # segment the methylation data
    # minSeg: Minimum number of CpGs in the segments
    # G: Number of segment clases to try
    res=methSeg(filtered.myobj[[1]],minSeg=10,G=1:6,
                join.neighbours = TRUE)

    plot(res$seg.mean,
     log10(width(res)),pch=20,
     col=scales::alpha(rainbow(length(unique(res$seg.group)))[as.numeric(res$seg.group)], 0.2),
     ylab="log10(length)",
     xlab="methylation proportion")

    # Annotation
    library(genomation)
    gene.obj=readTranscriptFeatures(system.file("extdata", "refseq.hg18.bed.txt", 
                                           package = "methylKit"))

    # creates a matrix containing percent methylation values
    perc.meth=percMethylation(meth)
}

Q1 <- factor(cells, levels=unique(cells))
as.numeric(Q1) - 1
c(1,2,3)[1:2]

