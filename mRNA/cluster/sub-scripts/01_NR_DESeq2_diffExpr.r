coldata <- 'Infer'
path <- '/home/jmendietaes/data/2021/mRNA/allProcessed/counts/salmon'
sample_name <- '/home/jmendietaes/data/2021/mRNA/allProcessed/furtherAnalysis/01_projectRestart/counts/tximportMerge'
cores <- 8
# Reference cell
refCell <- 'Mye'
#library(SummarizedExperiment)
library(tximport)
#library(tximeta)
library("DESeq2")
library(BiocParallel)


# args = commandArgs(trailingOnly=TRUE)
# if (length(args) < 2) {
#     stop("Usage: salmon_tximport.r <coldata> <salmon_out>", call.=FALSE)
# }

# coldata = args[1]
# path = args[2]
# sample_name = args[3]

prefix = sample_name
tx2gene = "salmon_tx2gene.tsv"
info = file.info(tx2gene)
if (info$size == 0) {
    tx2gene = NULL
} else {
    rowdata = read.csv(tx2gene, sep="\t", header = FALSE)
    colnames(rowdata) = c("tx", "gene_id", "gene_name")
    tx2gene = rowdata[,1:2]
}

fns = list.files(path, pattern = "quant.sf", recursive = T, full.names = T)
# I have the bad tendency to move things to folders called "old"
# This line is to skip those entries
fns <- Filter(function(x) !any(grepl("/old/", x)), fns)
names = basename(dirname(fns))
names(fns) = names
coldata = data.frame(files = fns, names = names)
## Get metadata
# fns <- dirname(fns)
# names = unlist(lapply(fns, function (x) strsplit(x, "/salmon/")[[1]][2]))
# names = unlist(lapply(names, function (x) strsplit(x, "/quant")[[1]][1]))
#names(fns) = names
coldata = data.frame(files = fns, names = names)
# Get batch and cel type metadata
samples.vec <- sort(rownames(coldata))
#groups <- sub("_.*$", "", samples.vec)
groups <- sub("^([^_]*_[^_]*).*", "\\1", samples.vec)
# First we get first element before "_" as cell type
cells <- sub("_.*", "\\1", groups)
# Second element after "_" is the batch
batches <- gsub("^.*_", "", groups)
# Add it to coldata
coldata["batch"] <- batches
coldata["cell"] <- cells


txi = tximport(fns, type = "salmon", tx2gene = tx2gene, 
                countsFromAbundance = "lengthScaledTPM")

# If only one batch
if (length(unique(coldata[,"batch"])) == 1) {
    dds <- DESeqDataSetFromTximport(txi, coldata, 
                            design = ~cell)
} else {
    print("Multiples batches not optimised yet")
    quit(save = "no", status = 0, runLast = FALSE)
    dds <- DESeqDataSetFromTximport(txi, coldata[,c("batch", "cell")], 
                            design = ~ batch + ko)
    normalized_counts <- counts(dds, normalized = TRUE)
    vsd <- vst(dds, blind = FALSE)
    mat <- assay(vsd)
    mat <- limma::removeBatchEffect(mat, vsd$batch)
    assay(vsd) <- mat
}

#-------------------------------------------
#- remove genes with <= 5 counts in all samples
#-------------------------------------------
# https://github.com/UCL-BLIC/DGE/blob/master/bin/run_deseq2.R
dds <- dds[ rowSums(counts(dds)) > 5, ]

# Set specific cell type as reference
dds$cell <- relevel(dds$cell, ref = refCell)


#---------
#- Run DGE
#----------
dds <- DESeq(dds, parallel=TRUE, BPPARAM=MulticoreParam(cores))
res <- results(dds)
summary(res)

pos <- res[,"padj"] <= 0.01
posFCup <- res[,"log2FoldChange"] >= 0.75
posFCdown <- res[,"log2FoldChange"] <= -0.75

sum(pos, na.rm=T)
sum(pos & posFCup, na.rm=T)
sum(pos & posFCdown, na.rm=T)

## Add gene name
idConv <- unique(rowdata[,c("gene_id", "gene_name")])
rownames(idConv) <- idConv[,"gene_id"]
res["gene_name"] <- idConv[rownames(res), "gene_name"]


if (!opt$vst) {
            rld <- rlog(dds)
        } else {
            rld <- vst(dds)
        }

        normalized_counts <- counts(dds, normalized = TRUE)
        vsd <- vst(dds, blind = FALSE)
        mat <- assay(vsd)
        mat <- limma::removeBatchEffect(mat, vsd$batch)
        assay(vsd) <- mat
        save(dds,rld,vsd, file=DDSFile)



# Create Tximeta object
se = SummarizedExperiment(assays = list(counts = txi[["counts"]], abundance = txi[["abundance"]], length = txi[["length"]]),
                        colData = DataFrame(coldata),
                        rowData = rowdata)
gi = summarizeToGene(txi, tx2gene = tx2gene)
gi.ls = summarizeToGene(txi, tx2gene = tx2gene,countsFromAbundance="lengthScaledTPM")


txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
# then below...
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)

## Preaprae for DESeq2

se <- tximeta(coldata)
gse <- summarizeToGene(se)
dds <- DESeqDataSet(gse, ~condition)

# Ensure that we have at least two cells
print(unique(cells))
if (length(unique(cells)) == 1) {
    quit(save = "no", status = 0, runLast = FALSE)
}

dds <- DESeqDataSetFromMatrix(countData = round(counts[, rownames(coldataSel) ]), 
                        colData = coldataSel[, c("ko", "batch")], design = ~ batch + cell)

# Load to DESeq2 
dds <- DESeqDataSet(gi.ls, design = ~ batch + cell)

# Using tximport as input
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ condition)


# Using Tximeta SummarizedExperiment as input
ddsTxi <- DESeqDataSet(se, design = ~ condition)
