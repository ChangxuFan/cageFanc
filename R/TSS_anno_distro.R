plotAnnoPie.fanc <- function(anno, out.file, main = NULL, ...) {
  system(paste0("mkdir -p ", dirname(out.file)))
  png(filename = out.file, width = 800, height = 500, res = 100, ...)
  try((ChIPseeker::plotAnnoPie(anno, main = main)))
  dev.off()
}

TSS.annotate.pipe <- function(df=NULL, subset.to, genome, plot.dir, just.ce = F,
                              promoter.boundary = 500, ce = NULL, debug =F,
                              genome.name = NULL, annoDb = NULL, TxDb = NULL) {
  # df should have 2 columns: bam, sample

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

  system(paste0("mkdir -p ", plot.dir))
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

bam2tss.fanc <- function(bam.df, genome=NULL, genome.name = NULL, out.dir, ce.only = F,
                         threads.master = 3, threads.sub = 3) {
  system(paste0("mkdir -p ", out.dir))
  if (is.character(bam.df))
    bam.df <- read.table(bam.df, sep = "\t", as.is = T)
  if (is.null(bam.df$genome))
    bam.df$genome <- NA

  if (is.null(bam.df$genome.name))
    bam.df$genome.name <- NA
  bam.df[bam.df == ""] <- NA
  ce.list <- bam.df %>% split(f=1:nrow(bam.df)) %>%
    mclapply(function(x) {
      if (is.null(genome))
        genome <- x$genome
      if (is.null(genome.name))
        genome.name <- x$genome.name

      if (is.na(genome.name) || is.null(genome.name)) {
        if (genome == "mm10") {
          genome.name <- "BSgenome.Mmusculus.UCSC.mm10"
        } else if (genome == "hg38") {
          genome.name <- "BSgenome.Hsapiens.UCSC.hg38"
        } else {
          stop("only mm10 and hg38 supported as arguments to genome parameter")
        }
      }
      ce <- CAGEr::CAGEexp( genomeName     = genome.name,
                            inputFiles     = x$bam,
                            inputFilesType = "bamPairedEnd",
                            sampleLabels   = x$sample)
      CAGEr::getCTSS(ce,correctSystematicG=FALSE,removeFirstG=FALSE,
                    mappingQualityThreshold = 0,
                     useMulticore = T, nrCores = threads.sub)
      # the additional G/C at the cap is regarded as softclipping by STAR, so don't
      #need to remove it.
      if (ce.only == T) {
        return(ce)
      }
      sample.out.dir <- paste0(out.dir, "/", x$sample, "/")
      system(paste0("mkdir -p ", sample.out.dir))
      trash <- bdg.bw.gen(object = ce, outdir = sample.out.dir, values = "raw")
      ucsc.bdg <- paste0(sample.out.dir, "/All.samples.CTSS.raw.bedGraph")
      cmd <- paste0("~/R_for_bash/ucscSplit.R -z ", " -i ", ucsc.bdg,
                    " -o ", sample.out.dir, " -t ", threads.sub)
      utilsFanc::cmd.exec.fanc(cmd = cmd, intern = F, run = T)

      return(ce)
    }, mc.cores = threads.master, mc.cleanup = T)

  names(ce.list) <- bam.df$sample
  saveRDS(ce.list, paste0(out.dir, "/ce.list.Rds"))
  return(ce.list)
}


peak.annotate.pipe <- function(ce, normalization = T,peak.call = T, peak.summary = T, plot =T,
                               genome, work.dir, plot.dir = NULL,
                               tpm.threshold = 2, nrPassThreshold = NULL, method = "distclu",
                               maxDist = 20, removeSingletons = T, keepSingletonsAbove = 4,
                               peak.watch.seed = 42, peak.watch.n = 100,
                               # merge.peak.min.samples = NULL,
                               peak.aggr.threashold = 5,
                               cap.filter.pct = 0.7,
                               promoter.boundary = 500, debug = F, save.Rds = NULL, ...) {

  if (is.null(plot.dir))
    plot.dir <- paste0(work.dir, "/plots")
  system(paste0("mkdir -p ", work.dir, " ", plot.dir))

  samples <- ce@colData$sampleLabels

  if (is.null(save.Rds))
    save.Rds <- paste0(work.dir, "/ce.Rds")

  # if (is.null(merge.peak.min.samples))
  #   merge.peak.min.samples <- length(samples)

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

  if (normalization == T) {
    png(filename = paste0(plot.dir, "/powerlaw.png"), width = 500, height = 500, res = 100)
    try(print(CAGEr::plotReverseCumulatives(ce, fitInRange = c(5, 1000), onePlot = TRUE)))
    dev.off()

    CAGEr::normalizeTagCount(ce, method = "powerLaw", fitInRange = c(5, 1000), alpha = 1.2, T = 10^6)
    # figure out what it means to have "1 tpm":
    counts.df <- assays(ce@ExperimentList$tagCountMatrix)$counts %>% as.data.frame()
    norm.df <- assays(ce@ExperimentList$tagCountMatrix)$normalizedTpmMatrix %>% as.data.frame()
    pos <- ce@ExperimentList@listData$tagCountMatrix %>% rowRanges() %>% as.data.frame() %>%
      mutate(pos = paste0(seqnames, ":", pos)) %>% pull(pos)
    # samples <- colnames(counts.df)
    count.norm.map <- lapply(samples, function(sample) {
      map.df <- data.frame(count =  counts.df[, sample, drop = T],
                           norm = norm.df[, sample, drop = T],
                           pos = pos) %>%
        filter(count > 0) %>%
        arrange(count, norm, pos)
      return(map.df)
    })
    names(count.norm.map) <- samples
    count.norm.map.2 <- lapply(samples, function(sample) {
      map.df <- data.frame(count =  counts.df[, sample, drop = T],
                           norm = norm.df[, sample, drop = T] %>% floor(),
                           pos = pos) %>%
        filter(count > 0) %>% group_by(norm) %>%
        summarise(!!(paste0("min_", sample)) := min(count)) %>% ungroup()
      return(map.df)
    })

    map.all.df <- Reduce(full_join, count.norm.map.2) %>% arrange(norm)

    saveRDS(count.norm.map, paste0(work.dir, "/count.norm.map.Rds"))
    saveRDS(map.all.df, paste0(work.dir, "/map.all.df.Rds"))
  }

  if (peak.call == T) {
    if (is.null(nrPassThreshold))
      nrPassThreshold <- samples %>% length()
    # browser()
    CAGEr::clusterCTSS( object = ce,
                        threshold = tpm.threshold,
                        thresholdIsTpm = T,
                        nrPassThreshold = nrPassThreshold,
                        method = method,
                        maxDist = maxDist,
                        removeSingletons = removeSingletons,
                        keepSingletonsAbove = keepSingletonsAbove, ...)

    saveRDS(ce, save.Rds)
  }

  if (peak.summary == T) {
    peaks.per.sample <- lapply(samples, function(sample) {
      df <- CAGEr::tagClusters(ce, sample) %>% as.data.frame()

      if (peak.watch.n > nrow(df))
        peak.watch.n <- nrow(df)
      set.seed(seed = peak.watch.seed)
      df.watch <- df[sample(1:nrow(df), size = peak.watch.n, replace = F),] %>% arrange(tpm)
      watch.file <- paste0(work.dir, "/tracks/", sample, "_peaks_watch.tsv")
      trash <- utilsFanc::write.zip.fanc(df = df.watch, out.file = watch.file, bed.shift = F, zip = F,
                                         col.names = T, row.names = F)

      df <- df %>%  mutate(forth = paste0(nr_ctss, ":", round(tpm), ":", round(tpm.dominant_ctss))) %>%
        dplyr::select(chr, start, end, forth, dominant_ctss, strand)
      track.file <- paste0(work.dir, "/tracks/", sample, "_peaks.bed")
      ##note: I turned off bed.shift. It seems that the tagcluster returned by cageR is already
      #bed format.
      trash <- utilsFanc::write.zip.fanc(df = df, out.file = track.file, bed.shift = F)

      peaks <- df %>% dplyr::select(chr, start, end)
      return(peaks)
    })
    # browser()

    sum.df <- data.frame(sample = samples,
                         n.peaks = sapply(peaks.per.sample, nrow) )
    print(sum.df)
    trash <- utilsFanc::write.zip.fanc(df = sum.df, out.file = paste0(work.dir, "/summary_n_peak_per_sample.tsv"),
                                       col.names = T, zip = F, row.names = F)

    #
    # names(peaks.per.sample) <- samples
    # reproducible.peaks <- lapply(samples, function(sample) {
    #   peaks.q <- peaks.per.sample[[sample]]
    #
    #
    #
    #   pass <- peaks.q %>% split(., f = 1:nrow(.)) %>%
    #     sapply(function(x) {
    #       overlap <- sapply(peaks.per.sample, function(peaks.s) {
    #         gr.q <- makeGRangesFromDataFrame(df = x)
    #         gr.s <- makeGRangesFromDataFrame(df = peaks.s)
    #         o <- findOverlaps(gr.q, gr.s)
    #         if (length(o) > 0)
    #           return(T)
    #         else
    #           return(F)
    #       })
    #       return(sum(overlap) >= merge.peak.min.samples)
    #     })
    #   names(pass) <- NULL
    #   return(peaks.q[pass, ])
    # })
    aggregateTagClusters(ce, tpmThreshold = peak.aggr.threashold, maxDist = maxDist)

    cc.df <- consensusClusters(ce, sample = NULL, returnInterquantileWidth = FALSE,
                               qLow = NULL, qUp = NULL) %>%
      mutate(consensus.cluster = NULL)
    cc.df <- cc.df %>% mutate(pct = rank(tpm)/nrow(cc.df)) %>%
      dplyr::select(chr, start, end, tpm, pct, strand)

    uG.pct <- cap.filter.nakul(df = cc.df, ce = ce, samples = samples, plot.dir = work.dir)
    bPass.cap <- uG.pct$uG.pct$max > cap.filter.pct
    # write(sum(bPass.cap), paste0(work.dir, "/nPassCapFilter_", cap.filter.pct, ".txt"))
    cc.df.cap <- cc.df[bPass.cap, ]
    # browser()
    metadata(ce)[[paste0("consensus_cap_", cap.filter.pct)]] <- cc.df.cap %>% makeGRangesFromDataFrame(keep.extra.columns = T)
    trash <- mapply(function(df, name) {
      trash <- sub.f.write.peaks(df = df, peak.watch.n = peak.watch.n,
                                 root.name = name,
                                 id = "tpm", seed = peak.watch.seed, bed.shift = T)
    }, list(cc.df, cc.df.cap), paste0(work.dir, "/tracks/consensus", c("", "_capfilter")))
    sum.df <- data.frame(type = c("pre-capfilter", "post-capfilter"),
                         n.peaks = c(nrow(cc.df), nrow(cc.df.cap)))
    write.table(sum.df, paste0(work.dir, "/summary_nPassCapFilter_", cap.filter.pct, ".tsv"),
                quote = F, col.names = T, row.names = F, sep = "\t")
    print("cap.filter.pct: " %>% paste0(cap.filter.pct))
    print(sum.df)
  }


  if (plot == T) {
    try(trash <- lapply(c("consensus", "consensus_capfilter", samples), function(x) {
      if (x == "consensus")
        gr <- makeGRangesFromDataFrame(cc.df)
      else if (x == "consensus_capfilter") {
        gr <- makeGRangesFromDataFrame(df = cc.df.cap)
      } else {
        gr <- CAGEr::tagClustersGR(ce)[[x]]
      }
      anno <- gr %>%
        ChIPseeker::annotatePeak(tssRegion=c(-1 * promoter.boundary, 0),
                                 TxDb=TxDb,
                                 annoDb=annoDb)
      trash <- plotAnnoPie.fanc(anno, out.file = paste0(plot.dir, "/", x, "_peaks_distro_pie.png"))
      return()
    }))
  }

  saveRDS(ce, save.Rds)
  return(ce)

}

sub.f.write.peaks <- function(df, peak.watch.n, root.name, id, seed, bed.shift) {
  track.file <- paste0(root.name, "_peaks.bed")
  ##note: the format returned by cager is very naughty.
  ## when you retrieve peaks for individual samples, they are already bed-shifted:
  #t <- tagClusters(ce, "B6_Ly49Dn_mRNA" )
  #min(t$end - t$start) # returns 1
  ## when you retrieve concensus peaks, they are not bed-shifted:
  #t <- consensusClusters(ce)
  #min(t$end - t$start) # this would give you 0
  trash <- utilsFanc::write.zip.fanc(df = df,
                                     out.file = track.file, bed.shift = bed.shift)
  if (peak.watch.n > nrow(df))
    peak.watch.n <- nrow(df)
  set.seed(seed =seed)
  df.watch <- df[sample(1:nrow(df), size = peak.watch.n, replace = F), ] %>%
    arrange(!!as.name(id))
  watch.file <- paste0(root.name, "_peaks_watch.tsv")
  trash <- utilsFanc::write.zip.fanc(df = df.watch,
                                     out.file = watch.file, bed.shift = bed.shift, zip = F,
                                     col.names = T)
  return()
}


filter.gr1.with.gr2 <- function(gr1, gr2) {
  if (is.data.frame(gr1))
    gr1 <- makeGRangesFromDataFrame(gr1, keep.extra.columns = T)
  if (is.data.frame(gr2))
    gr2 <- makeGRangesFromDataFrame(gr2, keep.extra.columns = T)

  o <- findOverlaps(query = gr1, subject = gr2)
  gr.out <- gr1[queryHits(o) %>% unique(),]
  return(gr.out)
}



gr.write.bed.fanc <- function(gr, out.file, use.strand = T, forth = NULL, fifth = NULL) {
  df <- as.data.frame(gr) %>% mutate(start = start - 1)
  if (!is.null(forth))
    df$forth <- df[, forth]
  else
    df$forth <- ""

  if (!is.null(fifth))
    df$fifth <- df[, fifth]
  else
    df$fifth <- ""

  if (use.strand == F)
    df$strand <- "."
  trash <- utilsFanc::write.zip.fanc(df = df %>% dplyr::select(seqnames, start, end, forth,fifth, strand),
                                     out.file = out.file, bed.shift = F, zip = T)
  return(out.file)
}

gr.read.bedgraph <- function(gr, bdg, format = "bedGraph", out.file = NULL) {
  if (is.data.frame(gr))
    gr <- makeGRangesFromDataFrame(gr)
  # browser()
  if (is.character(bdg)) {
    bdg <- rtracklayer::import(bdg, format = format, which = gr)
  }
  o <- findOverlaps(query = gr, subject = bdg)

  gr$queryHits <- 1:length(gr)
  bdg$subjectHits <- 1:length(bdg)
  value.df <- left_join(as.data.frame(o), as.data.frame(bdg)) %>% group_by(queryHits) %>%
    summarise(score = sum(score)) %>% ungroup() %>% as.data.frame()
  rm(bdg)
  rm(o)
  counts.df <- left_join(gr %>% `names<-`(NULL) %>% as.data.frame(), value.df)
  counts.df[is.na(counts.df)] <- 0
  if (!is.null(out.file)) {
    trash <- utilsFanc::write.zip.fanc(df = counts.df, out.file = out.file, bed.shift = T)
  }
  counts.gr <- makeGRangesFromDataFrame(counts.df, keep.extra.columns = T)
  return(counts.gr)
}

gr.read.bedgraph.m <- function(gr, bdg.files, format = "bedGraph",
                               sample.names = NULL, threads = 4,
                               factor = 1000000) {
  if (is.null(sample.names))
    sample.names <- names(bdg.files)
  if (is.null(sample.names))
    stop("samples must be named!!")
  names(bdg.files) <- sample.names
  counts.mat <- mclapply(sample.names, function(x) {
    file <- bdg.files[x]
    signal <- gr.read.bedgraph(gr = gr, bdg = file, format = format)$score
    return(signal)
  }, mc.cores = threads, mc.cleanup = T) %>% Reduce(cbind,.) %>%
   `colnames<-`(sample.names)
  counts.mat[is.na(counts.mat)] <- 0
  cpm.mat <- counts.mat %*% diag(1/colSums(counts.mat, na.rm = T)) * factor
  colnames(cpm.mat) <- colnames(counts.mat)
  mcols(gr) <- NULL
  se <- SummarizedExperiment(assays = list(counts = counts.mat, cpm = cpm.mat),
                             rowRanges = gr )
  return(se)
}
