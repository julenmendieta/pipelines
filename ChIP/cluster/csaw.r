library(csaw)
library(edgeR)


print(cbpdata)
## DataFrame with 4 rows and 3 columns
##          Name       Description      Path
##   <character>       <character>    <List>
## 1  SRR1145787 CBP wild-type (1) <BamFile>
## 2  SRR1145788 CBP wild-type (2) <BamFile>
## 3  SRR1145789 CBP knock-out (1) <BamFile>
## 4  SRR1145790 CBP knock-out (2) <BamFile>

# We set the minimum mapping quality score to 
# 10 to remove poorly or non-uniquely aligned reads
param <- readParam(minq=10, discard=blacklist)

## Computing the average fragment length
# For that we obtain the strand cross-correlaton 
x <- correlateReads(cbpdata$Path, param=reform(param, dedup=TRUE))
frag.len <- maximizeCcf(x)
plot(1:length(x)-1, x, xlab="Delay (bp)", ylab="CCF", type="l")
abline(v=frag.len, col="red")
text(x=frag.len, y=min(x), paste(frag.len, "bp"), pos=4, col="red")

## Count reads by windows
# In this case we use window size of 10bp
# Even if spacing is 50, the reads are longer than that, so reads 
#lying in the interval between adjacent windows will still be counted 
#into several windows
windowSize <- 10
win.data <- windowCounts(cbpdata$Path, param=param, 
                            width=windowSize,  # bin size
                            ext=frag.len,      # average length of sequence fragments
                            spacing=50,        # distance between consecutive windows
                            filter=10)         # min bin reads (from all bamfile) to exclude

## Filtering low abundance regions
filterBinSize <- 10000
bins <- windowCounts(cbpdata$Path, bin=TRUE, width=10000, param=param)
filter.stat <- filterWindowsGlobal(win.data, bins)

min.fc <- 3
keep <- filter.stat$filter > log2(min.fc)
print(paste0("Summary of filtered ", filterBinSize, "bins"))
summary(keep)
filtered.data <- win.data[keep,]

hist(filter.stat$filter, main="", breaks=50,
    xlab="Background abundance (log2-CPM)")
abline(v=log2(min.fc), col="red")


## Normalisation by TMM for composition biases

# ---- case of KO experiments ------
win.data <- normFactors(bins, se.out=win.data)
print('TMM normalisation factor scales:')
(normfacs <- win.data$norm.factors)

#I THINK THIS MIGHT ONLY BE FOR KO EXPERIMENTS, CHECK NORM FOR TRENDED BIASES

# Visualize effect of normalisation with mean-difference plots
# The dense cloud in each plot represents the majority of bins 
#in the genome. These are assumed to mostly contain background 
#regions. The red line represents the log-ratio of normalization 
#factors and passes through the centre of the cloud in each plot, 
#indicating that the bias has been successfully identified and removed

# NEED TO IMPLEMENT THE WAY TO STORE THIS PLOTS IN AN INFORME-LIKE
# AND DO ALL POSIBLE COMPARISONS
bin.ab <- scaledAverage(bins)
adjc <- calculateCPM(bins, use.norm.factors=FALSE)
print("Abundance-dependent trend in the log-fold change")
par(cex.lab=1.5, mfrow=c(1,3))
smoothScatter(bin.ab, adjc[,1]-adjc[,4], ylim=c(-6, 6),
    xlab="Average abundance", ylab="Log-ratio (1 vs 4)")
abline(h=log2(normfacs[1]/normfacs[4]), col="red")

smoothScatter(bin.ab, adjc[,2]-adjc[,4], ylim=c(-6, 6),
    xlab="Average abundance", ylab="Log-ratio (2 vs 4)")
abline(h=log2(normfacs[2]/normfacs[4]), col="red")

smoothScatter(bin.ab, adjc[,3]-adjc[,4], ylim=c(-6, 6),
    xlab="Average abundance", ylab="Log-ratio (3 vs 4)")
abline(h=log2(normfacs[3]/normfacs[4]), col="red")

# ---- case of diff cells and i hope peaks -----
# compute a matrix of offsets for model fitting
filtered.data <- normOffsets(filtered.data)
#head(assay(filtered.data, "offset"))

print("Abundance-dependent trend in the log-fold change after norm")
norm.adjc <- calculateCPM(filtered.data, use.offsets=TRUE)
norm.fc <- norm.adjc[,4]-norm.adjc[,1]
smoothScatter(win.ab, norm.fc, ylim=c(-6, 6), xlim=c(0, 5),
    xlab="Average abundance", ylab="Log-fold change")

lfit <- smooth.spline(norm.fc~win.ab, df=5)
lines(win.ab[o], fitted(lfit)[o], col="red", lty=2)

# The implicit assumption of non-linear methods is that most 
#windows at each abundance are not DB. Any systematic difference 
#between samples is attributed to bias and is removed. 
# The assumption of a non-DB majority is reasonable given that 
#the cell types being compared are quite closely related



## Statistical modelling
y <- asDGEList(filtered.data)
print("convert our RangedSummarizedExperiment object into a DGEList")
summary(y)

# Construct design matrix for experimental design
# TO BE CHANGED. in MY CASE WILL GO BY CHIP, IN THIS WAY WILL ALSO WORK
#WHEN I ADD THE REPLICATES
genotype <- cbpdata$Description
genotype[grep("wild-type", genotype)] <- "wt"
genotype[grep("knock-out", genotype)] <- "ko"

genotype <- factor(genotype)
design <- model.matrix(~0+genotype)
colnames(design) <- levels(genotype)
design

# estimate the negative binomial (NB) and quasi-likelihood (QL) 
#dispersions for each window 
y <- estimateDisp(y, design)
print("Modelling trend dispersion")
summary(y$trended.dispersion)
plotBCV(y)

print("Estimated prior d.f (Inf=all the QL dispersions are equal to the trend")
fit <- glmQLFit(y, design, robust=TRUE)
summary(fit$df.prior)
plotQLDisp(fit)

print("MDS plot with two dimensions for all samples")
plotMDS(cpm(y, log=TRUE), top=10000, labels=genotype,
    col=c("red", "blue")[as.integer(genotype)])


## Test for differential binding
# im going to have to check for new variables created that i will 
#combine pairwise as it is supposed to
contrast <- makeContrasts(wt-ko, levels=design)
res <- glmQLFTest(fit, contrast=contrast)

# Windows less than 100 bp apart are clustered into 
#regions with a maximum cluster width of finalResol bp
finalResol <- 5000
merged <- mergeResults(filtered.data, res$table, tol=100, 
    merge.args=list(max.width=finalResol))
print("Merged regions")
merged$regions
tabcom <- merged$combined
is.sig <- tabcom$FDR <= 0.05
print("Regions with FDR <= 0.05")
summary(is.sig)

print("Direction of change in significant bins")
table(tabcom$direction[is.sig])

# Direction according the best window in each cluster.
tabbest <- merged$best
is.sig.pos <- (tabbest$rep.logFC > 0)[is.sig]
summary(is.sig.pos)

# safe results to file
########### NEED TO ADD PATH FOR LOCATION
out.ranges <- merged$regions
mcols(out.ranges) <- DataFrame(tabcom,
    best.pos=mid(ranges(rowRanges(filtered.data[tabbest$rep.test]))),
    best.logFC=tabbest$rep.logFC)
saveRDS(file="cbp_results.rds", out.ranges)

