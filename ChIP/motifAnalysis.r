source(paste0(Sys.getenv("CODEBASE"), "sfb_jak_stat/src/00-init.R"))
require(LOLA)
require(doMC);registerDoMC(cores=SNOW.WORKERS)
ARGS <- commandArgs(trailingOnly=TRUE)

# Define parameters -------------------------------------------------
if(length(ARGS) < 3){
  # Manual settings
  peak.width <- NA # if null peaks / regions are not resized, else they are resized to the number provided
  bg.all.regions <- FALSE # if FALSE uses only significant regions as background, if TRUE uses all
  region.subtypes <- TRUE # if TRUE, splits region into promoter, intron, exon,..., if FALSE does not do that
} else {
  # From command line
  peak.width <- as.numeric(ARGS[1]) # if null peaks / regions are not resized, else they are resized to the number provided
  bg.all.regions <- ARGS[2] == TRUE # if FALSE uses only significant regions as background, if TRUE uses all
  region.subtypes <- ARGS[3] == TRUE # if TRUE, splits region into promoter, intron, exon,..., if FALSE does not do that
}
if(is.na(peak.width)) peak.width <- NULL

out.string <- paste(ifelse(is.null(peak.width), "NULL", peak.width), ifelse(bg.all.regions, "FullBG", "SigBG"), ifelse(region.subtypes, "SubtypeRegions", "AllRegions"), sep="_")
out.dir <- paste0("ATAC_12_Homer_", out.string, "/")
out <- dirout(paste0(out.dir))
flag(out)

out.lola <- dirout(paste0(out.dir, "/Lola/"))
out.homer <- dirout(paste0(out.dir, "/homer/"))
homer.map.file <- out("HomerMap.tsv")
jasparStatNames <- fread("metadata/JASPAR_Motifs.tsv")
jasparStatNames[, Name2 := paste0(Name, "-", ID)]

# FUNCTIONS ---------------------------------------------------------------
centerRegion <- function(gr, size){
  if(!is.null(size)){
    return(GenomicRanges::shift(resize(gr,size), width(gr)/2-size/2))
  } else {
    return(gr)
  }
}
grToBedFile <- function(gr){
  x <- data.table(
    chr=as.character(seqnames(gr)),
    start=start(gr),
    end=end(gr))
  x[, uniqueID := paste(chr, start, end, sep="_")]
  x$x <- 0
  x$Strand <- "*"
  x
}
combine.grlist <- function(grlist){x <- grlist[[1]]; for(i in 2:length(grlist)){x <- c(x, grlist[[i]])}; return(unique(x))}


# DATA --------------------------------------------------------------------
(load(PATHS$ANALYSIS$ATAC.DE))
(load(PATHS$CLEAN$ATAC$TPM))
ann <- PATHS$CLEAN$ATAC$Annotation


# MOTIFS ------------------------------------------------------------------
motifs.japsar.file <- dirout_Init("EXT_JASPAR_2_HOMER")("allMotifs.motif")
motifs.homer.file <- paste0(HOMER, "/data/knownTFs/vertebrates/known.motifs")
motifs.file <- out("all.motifs")
system(paste("cat", motifs.homer.file,motifs.japsar.file, ">", motifs.file))



# HOMER ANNOTATE REGIONS --------------------------------------------------
current.wd <- getwd()
setwd(out(""))

# Prepare regions
queryRegions <- toGR(ATAC.gmap$probe)
queryRegions <- centerRegion(queryRegions, size=peak.width)
allPeaks500 <- out("AllPeaks_Centered.bed")
write.tsv(grToBedFile(gr=queryRegions), file=allPeaks500, col.names=F)

# Run homer
if(!file.exists(homer.map.file)){
  system(paste(
    paste0(HOMER, "/bin/annotatePeaks.pl"),
    allPeaks500,
    "mm10",
    "-m", motifs.japsar.file,
    ">", homer.map.file))
}
setwd(current.wd)


# Define groups -----------------------------------------------------------

# Prepare annotations from HOMER
gmap.ann <- fread(homer.map.file)
colnames(gmap.ann)[1] <- "PeakID"
gmap.ann[, Annotation.clean := make.names(gsub(" \\(.+", "", Annotation))]

# Annotate differential peaks
ATAC.hits <- res
ATAC.hits[, direction := ifelse(logFC > 0, "up", "down")]
ATAC.hits <- merge(ATAC.hits, gmap.ann[,c("PeakID", "Annotation.clean"),with=F], by.x="rn", by.y="PeakID")
if(region.subtypes == FALSE) ATAC.hits$Annotation.clean <- "All"
ATAC.hits[,group := paste(cell_type, genotype, direction, make.names(Annotation.clean), sep="_")]

# Get significant peaks
ATAC.hits.sig <- ATAC.hits[padj < CUTOFFS$DE.ATAC$adj.P.Val]

# Plot significant peaks
ggplot(ATAC.hits.sig, aes(x=genotype, fill=factor(sign(logFC)))) +
  geom_bar(position="dodge") + xRot() + facet_grid(cell_type ~ .) +
  scale_y_log10()
ggsave(out("NumberOfHits.pdf"), w=6,h=12)

# Plot by region type
p <- ggplot(ATAC.hits.sig, aes(x=Annotation.clean)) + geom_bar() +
  facet_grid(cell_type + direction ~ genotype) +
  scale_y_log10() + theme_bw(12) + xRot()
ggsave(out("NumberOfHits_byRegionType.pdf"), w=15,h=24, plot=p)

p <- ggplot(ATAC.hits.sig, aes(x=cell_type)) + geom_bar() +
  facet_grid(Annotation.clean ~  genotype +  direction) +
  scale_y_log10() + theme_bw(12) + xRot()
ggsave(out("NumberOfHits_byRegionType2.pdf"), w=24,h=15, plot=p)


ATAC.hits.sig <- ATAC.hits.sig[order(abs(logFC), decreasing=T)][, head(.SD, n=500), by="group"]

# Plot by region type
p <- ggplot(ATAC.hits.sig, aes(x=Annotation.clean)) + geom_bar() +
  facet_grid(cell_type + direction ~ genotype) +
  scale_y_log10() + theme_bw(12) + xRot()
ggsave(out("NumberOfHits_byRegionType_AfterFiltering.pdf"), w=15,h=24, plot=p)

# Define background
bg.data <- if(bg.all.regions) ATAC.hits else ATAC.hits.sig
bg.by.ct.string <- lapply(with(bg.data, split(rn, paste(cell_type, Annotation.clean, sep="_"))), unique)
bg.by.ct <- lapply(bg.by.ct.string, toGR)

sort(sapply(bg.by.ct.string, length))


# HOMER -------------------------------------------------------------------
homer.file <- out.homer("homer.tsv")
if(!file.exists(homer.file)){
  nam <- ATAC.hits.sig$group[1]
  foreach(nam = unique(ATAC.hits.sig$group)) %dopar% {
    message(nam)
    if(!file.exists(out.homer(nam, "/", "knownResults.txt"))){
      dir.create(out.homer(nam))
      
      bg.grp <- paste(gsub("_.+$", "", nam), gsub(".+_(.+)$", "\\1", nam), sep="_")
      bg.regions <- bg.by.ct.string[[bg.grp]]
      hit.regions <- ATAC.hits.sig[group == nam]$rn
      stopifnot(all(hit.regions %in% bg.regions))
      
      homer.bg <- out.homer(nam, "/bg.bed")
      write.table(grToBedFile(toGR(bg.regions)), file=homer.bg, sep="\t", quote=F, row.names=F,col.names=F)
      
      homer.input <- out.homer(nam, "/",nam, ".bed")
      write.table(grToBedFile(toGR(hit.regions)), file=homer.input, sep="\t", quote=F, row.names=F,col.names=F)
      
      peak.width.arg <- if(is.null(peak.width)) "-size given" else paste("-size",peak.width)
      
      system(paste(
        paste0(HOMER, "/bin/findMotifsGenome.pl"),
        homer.input,
        "mm10",
        out.homer(nam, "/"),
        peak.width.arg,
        "-nomotif",
        "-mknown", motifs.file,
        "-bg", homer.bg,
        "&>", out.homer(nam, "/", "homer.log")))
    }
  }
  
  homerRes <- data.table()
  for(nam in unique(ATAC.hits.sig$group)){
    fx <- out.homer(nam, "/", "knownResults.txt")
    if(!file.exists(fx)) next
    x <- fread(fx, check.names=T)
    #colnames(x) <- make.names(gsub("\\(of \\d+\\)", "_count", colnames(x)))
    x[,targetPercent := as.numeric(gsub("%", "", X..of.Target.Sequences.with.Motif))]
    x[,bgPercent := as.numeric(gsub("%", "", X..of.Background.Sequences.with.Motif))]
    x[, log2Odds := log2(targetPercent/bgPercent + 1)]
    x$group <- nam
    homerRes <- rbind(homerRes, x, fill=TRUE)
  }
  homerRes$qval <- homerRes$q.value..Benjamini.
  homerRes$Gene <- gsub("^(.)", "\\U\\1", tolower(gsub("(.+?)[\\/|\\(].+", "\\1", homerRes$Motif.Name)), perl=TRUE)
  i <- 1
  for(i in 1:nrow(jasparStatNames)){
    homerRes[toupper(Gene) == jasparStatNames[i]$ID, Gene := jasparStatNames[i]$Name2]
  }
  x <- unique(homerRes[,c("Gene", "Motif.Name")])
  x[,Gene2 := make.unique(Gene)]
  homerRes$Gene2 <- NULL
  homerRes <- merge(homerRes, x, by=c("Gene", "Motif.Name"))
  write.tsv(homerRes, file=homer.file)
} else {
  homerRes <- fread(homer.file)
}

# PLOT
# remove duplicated motifs
duplicated.motifs <- unique(homerRes[, .N, by=c("Motif.Name", "group")][N > 1]$Motif.Name)
homerRes <- homerRes[!Motif.Name %in% duplicated.motifs]
homerRes[,Motif.short := substr(homerRes$Motif.Name,0,10)]
homerRes$Motif <- paste0(substr(homerRes$Motif.Name,0,10), "_", abbreviate(homerRes$Motif.Name, minlength=10), "_", homerRes$Consensus)
homerRes <- homerRes[, -grep("X..of.Background.Sequences.with.Motif", colnames(homerRes)), with=F]
homerRes <- homerRes[, -grep("X..of.Target.Sequences.with.Motif.of", colnames(homerRes)), with=F]

pDT.all <- copy(homerRes)
pDT.all <- cbind(
  pDT.all,
  setNames(data.table(do.call(rbind, strsplit(pDT.all$group, "_"))), c("cell_type", "genotype", "direction", "region_type"))
)
pDT.all[log2Odds == Inf, log2Odds := max(pDT.all[log2Odds != Inf]$log2Odds)]
pDT <- copy(pDT.all)
str(freqMotifs <- unique(as.character(pDT[qval < 0.05][order(log2Odds)]$Gene2))[1:100])
pDT <- hierarch.ordering(pDT, toOrder="Gene2", orderBy="group", value.var="log2Odds")
pDT.homer.all <- pDT[P.value < 0.9][Gene2 %in% freqMotifs]
if(nrow(pDT.homer.all) > 2){
  p <- ggplot(pDT.homer.all,
              aes(x=region_type, y=Gene2, color=log2Odds, size=pmin(5, -log10(qval)))) +
    geom_point() + theme_bw(12) +
    facet_grid(cell_type + direction ~ genotype, scales="free", space="free") +
    xRot() +
    scale_color_gradient2(low="blue", mid="black", high="red")
  ggsave(out("Homer_HM_top100",".pdf"), width=25, height=120, limitsize=F, plot=p)
  
  p <- ggplot(pDT.homer.all[, .(N = length(unique(Gene2))), by=c("region_type", "cell_type", "direction", "genotype")],
              aes(x=region_type, y=N)) +
    facet_grid(cell_type + direction ~ genotype) +
    theme_bw(12) + xRot() +
    geom_bar(stat="identity")
  ggsave(out("Homer_HM_top100_Counts.pdf"), w=15,h=24, plot=p)
}

pDTxxx <- pDT[Gene2 %in% jasparStatNames$Name2]
if(nrow(pDTxxx) > 0){
  p <- ggplot(pDT[Gene2 %in% jasparStatNames$Name2],
              aes(x=region_type, y=Gene2, color=log2Odds, size=pmin(5, -log10(qval)))) +
    geom_point() + theme_bw(12) +
    facet_grid(cell_type + direction ~ genotype, scales="free", space="free") +
    xRot() +
    scale_color_gradient2(low="blue", mid="black", high="red")
  ggsave(out("Homer_HM_STAT",".pdf"), width=25, height=30, limitsize=F, plot=p)
}



# LOLA --------------------------------------------------------------------
lola.file <- out.lola("lola.tsv")
ATAC.hits.sig[,lola.grp := paste(cell_type, Annotation.clean, sep="_")]
if(!file.exists(lola.file)){
  lola.res <- data.table()
  LOLA_regionDB = loadRegionDB(dbLocation=LOLA.CORE)
  lola.grpx <- "M_intron"
  for(lola.grpx in unique(ATAC.hits.sig$lola.grp)){
    lola.ct.hits <- lapply(with(ATAC.hits.sig[lola.grp == lola.grpx], split(rn, group)), toGR)
    lola.res.ct <- runLOLA(userSets=lapply(lola.ct.hits, function(x) centerRegion(x, peak.width)), userUniverse=centerRegion(bg.by.ct[[lola.grpx]], peak.width), regionDB=LOLA_regionDB)
    lola.res.ct <- cleanLOLA(lola.res.ct)
    lola.res <- rbind(lola.res, data.table(lola.res.ct))
  }
  x <- unique(lola.res[,c("target", "filename")])
  x[,target2 := make.unique(target)]
  lola.res <- merge(lola.res, x, by=c("target", "filename"))
  lola.res[,logOddsRatio := log2(oddsRatio)]
  lola.res[logOddsRatio == -Inf, logOddsRatio := min(lola.res[logOddsRatio != -Inf]$logOddsRatio)]
  lola.res[logOddsRatio == Inf, logOddsRatio := max(lola.res[logOddsRatio != Inf]$logOddsRatio)]
  lola.res <- cbind(
    lola.res,
    setNames(data.table(do.call(rbind, strsplit(as.character(lola.res$userSet), "_"))), c("cell_type", "genotype", "direction", "region_type"))
  )
  write.tsv(lola.res, lola.file)
} else {
  lola.res <- fread(lola.file)
}

lola.sig <- lola.res[BH < 0.1]
lola.split.factors <- with(lola.sig[!is.na(target)], split(target, factor(userSet)))
if(length(lola.split.factors) > 0){
  all.factors <- unique(do.call(c, lola.split.factors))
  gMat <- sapply(lola.split.factors, function(x) all.factors %in% x) + 0
  row.names(gMat) <- all.factors
  
  interesting.targets <- c(names(sort(colSums(t(gMat) * (1-colSums(gMat)/nrow(gMat))), decreasing=T)[1:35]), "GATA1", "GATA2")
  
  cleanDev()
  pdf(out("LOLA_FactorSummary.pdf"), width=12, height=12, onefile=F)
  pheatmap(gMat, fontsize_row=4)
  dev.off()
  
  cleanDev()
  pdf(out("LOLA_FactorSummary_Interesting.pdf"), width=12, height=12, onefile=F)
  pheatmap(gMat[interesting.targets,sort(colnames(gMat))])
  dev.off()
  
  xxx <- "hits"
  for(xxx in c("all", "hits")){
    lola.res.sig <- lola.res[BH < 0.1]
    if(xxx == "hits") lola.res.sig <- lola.res.sig[target %in% interesting.targets]
    if(nrow(lola.res.sig) < 2) next
    lola.height <- if(xxx=="hits") 50 else 100
    #if(xxx == "chemokines") lola.res.sig <- lola.res.sig[target %in% chemokine.genes]
    lola.res.sig <- lola.res.sig[!(is.na(target) | target == "NA")]
    lola.res.sig <- lola.res.sig[, .(logOddsRatio = max(logOddsRatio), pValueLog = max(pValueLog)), by=c("userSet", "target")]
    lola.res.sig <- cbind(
      lola.res.sig,
      setNames(data.table(do.call(rbind, strsplit(as.character(lola.res.sig$userSet), "_"))), c("cell_type", "genotype", "direction", "region_type"))
    )
    ggplot(lola.res.sig, aes(x=region_type, y=target, color=pmin(5, logOddsRatio), size=pmin(10, pValueLog))) +
      geom_point() + theme_bw(12) +
      #theme(axis.text.y = element_text(size=6)) +
      xRot() +
      scale_color_gradient2(low="blue", mid="blue", high="red") +
      scale_size_continuous(range=c(0.5,3)) +
      facet_grid(cell_type + direction ~ genotype, space="free", scales="free")
    #theme(strip.text.y=element_text(angle=0))
    ggsave(out("LOLA_RES_",xxx,".pdf"), w=14, h=50, limitsize=F)
  }
}



flag.rm(out)