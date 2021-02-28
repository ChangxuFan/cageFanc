plotAnnoPie.fanc <- function(anno, out.file, ...) {
  png(filename = out.file, width = 800, height = 500, res = 100, ...)
  try(print(ChIPseeker::plotAnnoPie(anno)))
  dev.off()
}

TSS.annotate.pipe <- function(df=NULL, subset.to, genome, plot.dir, just.ce = F,
                              promoter.boundary = 500, ce = NULL, debug =F) {
  # df should have 2 columns: bam, sample
  system(paste0("mkdir -p ", plot.dir))

  if (genome == "mm10") {
    genome.name <- "BSgenome.Mmusculus.UCSC.mm10"
    annoDb <- "org.Mm.eg.db"
    TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
  }

  if (genome == "hg38") {
    genome.name <- "BSgenome.Hsapiens.UCSC.hg38"
    annoDb <- "org.Hs.eg.db"
    TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
  }


  if (is.null(ce)) {
    if (is.character(df))
      df <- read.table(df, as.is = T, sep = "\t", header = T)
    ce <- CAGEr::CAGEexp( genomeName     = genome.name,
                   inputFiles     = df$bam,
                   inputFilesType = "bamPairedEnd",
                   sampleLabels   = df$sample)
    CAGEr::getCTSS(ce,correctSystematicG=FALSE,removeFirstG=FALSE,
            useMulticore = T, nrCores = nrow(df))
  }

  if (just.ce == T) {
    return(ce)
  }


  threads <- ce@sampleMap$primary %>% length()
  if (debug == T)
    threads <- 1
  mclapply(ce@sampleMap$primary, function(x) {
    gr <- CAGEr::CTSStagCountGR(ce, samples = x)
    chosen <- sample(1:length(gr), size = subset.to, replace = T, prob = gr$score %>%  as.numeric( ) )
    gr <- gr[chosen,]
    # here replace = T is a desired property. This command would give the effect of subsampling the reads
    #instead of the regions.
    gr <- GenomicRanges::sort(gr)
    anno <- ChIPseeker::annotatePeak(gr, tssRegion=c(-1 * promoter.boundary, 0),
                 TxDb=TxDb,
                 annoDb=annoDb)
    # Browse[2]> t <- anno@anno %>% as.data.frame()
    # Browse[2]> t$chrpos <- paste0(t$seqnames, ":", t$pos)
    # Browse[2]> t$chrpos %>% duplicated() %>% which() %>% length()
    # [1] 10455 # this means that the anno pipeline respects the fact that one TSS can have multiple reads
    trash <- plotAnnoPie.fanc(anno = anno, out.file = paste0(plot.dir, "/", x, "_reads_distro_pie.png"))
    return()
  }, mc.cores = threads)
  return(ce)
}




peak.annotate.pipe <- function(ce, skip.peak.call=F, genome, plot.dir, skip.plot =F,
                               skip.normalization = F,
                               promoter.boundary = 500, debug = F, save.Rds = NULL) {
  system(paste0("mkdir -p ", plot.dir))

  if (genome == "mm10") {
    genome.name <- "BSgenome.Mmusculus.UCSC.mm10"
    annoDb <- "org.Mm.eg.db"
    TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
  }

  if (genome == "hg38") {
    genome.name <- "BSgenome.Hsapiens.UCSC.hg38"
    annoDb <- "org.Hs.eg.db"
    TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
  }

  if (skip.normalization == F) {
    png(filename = paste0(plot.dir, "/powerlaw.png"), width = 500, height = 500, res = 100)
    try(print(CAGEr::plotReverseCumulatives(ce, fitInRange = c(5, 1000), onePlot = TRUE)))
    dev.off()

    CAGEr::normalizeTagCount(ce, method = "powerLaw", fitInRange = c(5, 1000), alpha = 1.2, T = 10^6)
  }

  if (skip.peak.call==F) {

    CAGEr::clusterCTSS( object = ce
                        , threshold = 1
                        , thresholdIsTpm = TRUE
                        , nrPassThreshold = 1
                        , method = "distclu"
                        , maxDist = 20
                        , removeSingletons = TRUE)
    if (!is.null(save.Rds))
      saveRDS(ce, save.Rds)
  }

  samples <- CAGEr::tagClustersGR(ce.fanc) %>% names()

  threads <- samples %>% length()
  if (debug == T)
    threads <- 1

  if (skip.plot == F) {
    try(trash <- mclapply(samples, function(x) {
      anno <- CAGEr::tagClustersGR(ce)[[x]] %>%
        ChIPseeker::annotatePeak(tssRegion=c(-1 * promoter.boundary, 0),
                                 TxDb=TxDb,
                                 annoDb=annoDb)
      plotAnnoPie.fanc(anno, out.file = paste0(plot.dir, "/", x, "_peaks_distro_pie.png"))
      return()
    }, mc.cores = threads, mc.cleanup = T))
  }

  return(ce)

}

# annotate.pipe.list <- function(df=NULL, subset.to, genome, work.dir, plot.dir,
#                                ce.list = NULL, plot.TSS.enrich = T, normalize = T,
#                                call.peaks = T,
#                                promoter.boundary = 500, debug =F, save.rds = T) {
#   system(paste0("mkdir -p ", plot.dir))

#   if (genome == "mm10") {
#     genome.name <- "BSgenome.Mmusculus.UCSC.mm10"
#     annoDb <- "org.Mm.eg.db"
#     TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
#   }

#   if (genome == "hg38") {
#     genome.name <- "BSgenome.Hsapiens.UCSC.hg38"
#     annoDb <- "org.Hs.eg.db"
#     TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
#   }

#   if (is.null(ce.list)) {
#     if (is.character(df))
#       df <- read.table(df, as.is = T, sep = "\t", header = T)
#     threads <- nrow(df)
#     if (debug == T)
#       threads <- 1
#     ce.list <- df %>% split(., f = factor(.$sample, levels = unique(.$sample))) %>%
#       mclapply(function(x) {
#         ce <- CAGEr::CAGEexp( genomeName     = genome.name,
#                         inputFiles     = x$bam,
#                         inputFilesType = "bamPairedEnd",
#                         sampleLabels   = x$sample)
#         CAGEr::getCTSS(ce,correctSystematicG=FALSE,removeFirstG=FALSE,
#                        useMulticore = F)
#         return(ce)
#       }, mc.cores = threads, mc.cleanup = T)
#     if (save.rds == T)
#       saveRDS(ce.list, paste0(work.dir, "/ce_list_1.Rds"))
#   }

#   threads <- length(ce.list)
#   if (debug == T)
#     threads <- 1



#   if (normalize == T) {
#     ce.list <- mapply(function(ce, name) {
#       png(filename = paste0(plot.dir, "/", name,"_powerlaw.png"), width = 500, height = 500, res = 100)
#       try(print(CAGEr::plotReverseCumulatives(ce, fitInRange = c(5, 1000), onePlot = TRUE)))
#       dev.off()

#       CAGEr::normalizeTagCount(ce, method = "powerLaw", fitInRange = c(5, 1000), alpha = 1.2, T = 10^6)
#       return(ce)
#     }, ce.list, names(ce.list), SIMPLIFY = F)

#   }

#   if (call.peaks ==T) {

#     ce.list <- mclapply

#     CAGEr::clusterCTSS( object = ce
#                         , threshold = 1
#                         , thresholdIsTpm = TRUE
#                         , nrPassThreshold = 1
#                         , method = "distclu"
#                         , maxDist = 20
#                         , removeSingletons = TRUE)
#     if (!is.null(save.Rds))
#       saveRDS(ce, save.Rds)
#   }

#   samples <- CAGEr::tagClustersGR(ce.fanc) %>% names()


#   }


# }
