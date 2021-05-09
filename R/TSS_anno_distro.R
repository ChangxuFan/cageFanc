plotAnnoPie.fanc <- function(anno, out.file, ...) {
  png(filename = out.file, width = 800, height = 500, res = 100, ...)
  try(print(ChIPseeker::plotAnnoPie(anno)))
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
                     useMulticore = T, nrCores = threads.sub)
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
    metadata(ce)[[paste0("consensus_cap_", cap.filter.pct)]] <- cc.df.cap
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
   `colnames<-`(names(atac.insert.sites))
  counts.mat[is.na(counts.mat)] <- 0
  cpm.mat <- counts.mat %*% diag(1/colSums(counts.mat, na.rm = T)) * factor
  colnames(cpm.mat) <- colnames(counts.mat)
  mcols(gr) <- NULL
  se <- SummarizedExperiment(assays = list(counts = counts.mat, cpm = cpm.mat), 
                             rowRanges = gr )
  return(se)
}

range.cor <- function(gr.q = NULL, gr.s = NULL, cage.gr=NULL, make.consensus = F,
                      cage.samples, atac.se = NULL,
                      atac.samples = NULL,
                      col.x, col.y, add.nearest.gene = F, rm.blacklist = F,
                      genome, annoDb, TxDb,  ...) {
  
  if (!is.null(cage.gr)) {
    if (make.consensus == T) {
      tpm.df <- mcols(cage.gr)[, cage.samples] %>% as.data.frame()
      tpm.df[is.na(tpm.df)] <- 0
      cage.gr$tpm <- tpm.df %>% 
        apply(1, mean)
    }
    gr.q <- cage.gr
  }
    
  if (!is.null(atac.se)) {
    if (is.null(atac.samples))
      atac.samples <- colnames(atac.se)
    gr.s <- rowRanges(atac.se)
    mcols(gr.s) <- data.frame(atac.cpm = assays(atac.se)$cpm[, atac.samples] %>% apply(1, mean))
  }
  # cabrowser()
  o.gr <- plyranges::join_overlap_left(x = gr.q, y = gr.s)  
  names(o.gr) <- NULL
  if (add.nearest.gene == T) {
    if (genome == "mm10") {
      genome.name <- "BSgenome.Mmusculus.UCSC.mm10"
      annoDb <- "org.Mm.eg.db"
      #annoDb <- NULL
      # TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
      TxDb <- loadDb(file = "~/genomes/mm10/gencode/TxDb.mm10.gencode.v24.klraps.sqlite")
    }
    
    if (genome == "hg38") {
      genome.name <- "BSgenome.Hsapiens.UCSC.hg38"
      # annoDb <- "org.Hs.eg.db"
      TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    }
    # browser()
    anno <- o.gr %>%
      ChIPseeker::annotatePeak(TxDb=TxDb, annoDb = annoDb)
    o.gr <- anno@anno
  }

  if (rm.blacklist  == T) {
    bl <- rtracklayer::import.bed(paste0("~/genomes/", genome, "/blacklist/", genome, ".blacklist.bed"))
    o.bl <- findOverlaps(o.gr, bl)
    o.gr <- o.gr[!(1:nrow(o.gr)) %in% queryHits(o.bl),]
  }
  
  o.df <- o.gr %>% as.data.frame()
  
  o.df[, col.x][is.na(o.df[, col.x])] <- 0
  o.df[, col.y][is.na(o.df[, col.y])] <- 0
  
  o.df <- o.df %>% mutate(locus = paste0(seqnames, ":", start, "-", end)) %>% 
    mutate(plotly = paste0(locus, "_", SYMBOL)) %>% 
    mutate(SYMBOL.ez = gsub("[a-zA-Z]","",SYMBOL) %>% sub("\\-", "",.))
  o.df$SYMBOL.ez[!grepl("Klra\\d",o.df$SYMBOL)] <- NA
  p <- scFanc::xy.plot(df = o.df, x = col.x, y = col.y,  ...)
  
  res <- list(p = p, o.gr = o.gr, o.df = o.df)
  return(res)
}

sync.atac.2.cage <- function(cage.gr, atac.se, atac.signal.files, genome, sample.names = NULL,
                             atac.signal.format, threads = 4, factor = 1000000,
                             region.size = 750, rm.black.list = T ) {
  atac.called.gr <- rowRanges(atac.se)
  mcols(atac.called.gr)$atac.called <- 1
  o.gr <- plyranges::join_overlap_left(x = cage.gr, y = atac.called.gr)  
  if (rm.black.list == T) {
    bl <- rtracklayer::import.bed(paste0("~/genomes/", genome, "/blacklist/", genome, ".blacklist.bed"))
    o.bl <- findOverlaps(o.gr, bl)
    o.gr <- o.gr[!(1:length(o.gr)) %in% queryHits(o.bl),]
    
    c.bl <- findOverlaps(cage.gr, bl)
    cage.gr <- cage.gr[!(1:length(cage.gr)) %in% queryHits(c.bl),]
  }
  re.read.id <- which(is.na(o.gr$atac.called))
  gr.re.read <- o.gr[re.read.id , ] %>% as.data.frame() %>% 
    mutate(mid = floor((start + end)/2)) %>% 
    mutate(start = floor(mid - region.size/2), end = floor(mid + region.size/2)) %>% 
    makeGRangesFromDataFrame()
  if (length(gr.re.read) == 0) {
    return(list(cage.gr = cage.gr, atac.se = atac.se))
  }
  if (is.null(sample.names))
    sample.names <- names(atac.signal.files)
  if (is.null(sample.names))
    stop("samples must be named!!")
  
  samples.ori <- colnames(atac.se)
  if (!identical(sample.names, samples.ori)) {
    stop("sample.names must exactly match the colnames in atac.se. In the same order")
  }
  
  re.read.se <- gr.read.bedgraph.m(gr = gr.re.read, bdg.files = atac.signal.files,
                                   sample.names = sample.names, threads = threads,
                                   format = atac.signal.format)
  new.atac.se <- rbind(atac.se, re.read.se)
  new.atac.se <- se.fast.normalize(se = new.atac.se, factor = factor, slot.name = "cpm")
  return(list(cage.gr = cage.gr, atac.se = new.atac.se))
}

se.fast.normalize <- function(se, factor = 1000000, slot.name = "cpm") {
  counts.mat <- assays(se)$counts
  assays(se)[[slot.name]] <- counts.mat %*% diag(1/colSums(counts.mat, na.rm = T)) * factor
  return(se)
}

cager.get.counts <- function(ce, gr, consensus = T, samples = NULL) {
  if (is.null(samples))
    samples <- ce$sampleLabels
  else
    samples <- ce$sampleLabels %>% .[.%in% samples]
  mcols(gr) <- NULL
  names(gr) <- NULL
  tpm.gr <- rowRanges(ce@ExperimentList$tagCountMatrix)
  tpm.mat <- assays(ce@ExperimentList$tagCountMatrix)$normalizedTpmMatrix 
  mcols(tpm.gr) <- tpm.mat
  gr$id <- 1:length(gr)
  j <- plyranges::join_overlap_left(gr, tpm.gr) %>% as.data.frame() %>% 
    group_by(id) %>% summarise_at(.vars = samples, .funs = sum) %>% as.data.frame()
  if (consensus == T) {
    j$tpm <- rowSums(j %>% mutate(id = NULL))
    j <- j %>% select(id, tpm)
  } 
  mcols(gr) <- left_join(mcols(gr) %>% as.data.frame(), j) %>% mutate(id = NULL)
  return(gr)

}

sync.atac.2.cage.2 <- function(cage.se, atac.se, genome, annotate = T, rmbl = T, 
                               bw.sample.map, bw.ext.size = 500, atac.signal.format = "BigWig",
                               factor = 100000,
                               threads) {
  ###
  #note: the single most important parameter is the id for j
  #it must be equal to the row numbers of j
  ### too many things depending on this!!
  cage.rowR <- rowRanges(cage.se) %>% plyranges::mutate(., cage.id = 1:length(.))
  atac.rowR <- rowRanges(atac.se) %>% plyranges::mutate(., atac.id = 1:length(.))
  j <- plyranges::join_overlap_left(x = cage.rowR,
                                    y = atac.rowR)
  
  if (annotate == T) {
    if (genome == "mm10") {
      genome.name <- "BSgenome.Mmusculus.UCSC.mm10"
      annoDb <- "org.Mm.eg.db"
      #annoDb <- NULL
      # TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
      TxDb <- AnnotationDbi::loadDb(file = "~/genomes/mm10/gencode/TxDb.mm10.gencode.v24.klraps.sqlite")
    }
    
    if (genome == "hg38") {
      genome.name <- "BSgenome.Hsapiens.UCSC.hg38"
      # annoDb <- "org.Hs.eg.db"
      TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    }
    # browser()
    anno <- j %>%
      ChIPseeker::annotatePeak(TxDb=TxDb, annoDb = annoDb)
    if (length(anno@anno) != length(j)) {
      stop ("length(anno@anno) != length(j)")
    }
    j <- anno@anno
    rm(anno)
  }
  
  if (rmbl == T) {
    bl <- rtracklayer::import.bed(paste0("~/genomes/", genome, "/blacklist/", genome, ".blacklist.bed"))
    j <- plyranges::filter_by_non_overlaps(x = j, y = bl)
  }
  
  j$id <- 1:length(j)

  re.read.id <- which(is.na(j$atac.id))
  gr.re.read <- j[re.read.id , ] + bw.ext.size
  
  if (length(gr.re.read) > 0) {
    samples.ori <- colnames(atac.se)
    bw.map <- read.table(bw.sample.map, header = T)
    # attempt to rearrange:
    bw.map <- bw.map[match(samples.ori, bw.map$sample),]
    
    if (!identical(bw.map$sample, samples.ori)) {
      stop("sample.names must exactly match the colnames in atac.se. In the same order")
    }
    
    re.read.se <- gr.read.bedgraph.m(gr = gr.re.read, bdg.files = bw.map$bw,
                                     sample.names = bw.map$sample, threads = threads,
                                     format = atac.signal.format)
    
    if (!identical(ranges(gr.re.read), ranges(rowRanges(re.read.se)))) {
      stop("identical(ranges(gr.re.read), ranges(rowRanges(re.read.se))) != T")
    }
    
    j$atac.id[re.read.id] <- (1:nrow(re.read.se)) + nrow(atac.se)
    atac.se <- rbind(atac.se, re.read.se)
    # atac.se <- se.fast.normalize(se = atac.se, factor = factor, slot.name = "cpm")
  }
  
  sync.o <- list(j = j, cage = cage.se[j$cage.id, ],
                 atac = atac.se[j$atac.id, ])
  
  
  sync.o$atac <- se.fast.normalize(se = sync.o$atac, factor = factor, slot.name = "cpm")
  
  return(sync.o)
}




sync.o.rescue <- function (sync.o = sync.o, factor = 100000,
                           bRescue, regex.rescue = NULL, rescue.col.name,
                           bam.sample.map, 
                           bam.ext.size,  rescue.shift,
                           tag, header.file, m.score.cutoff,
                           threads = NULL, rescue.report.dir) {
  ###
  #note: the single most important parameter is the id for j
  #it must be equal to the row numbers of j
  ### too many things depending on this!!
  ###
  ###
  j <- sync.o$j
  if (!is.null(regex.rescue)) {
    # note this is corresponing to the ranges of original cage!
    bRescue <- mcols(j)[, rescue.col.name] %>% grepl(regex.rescue, .)
  }
  # make it correspond to j:
  if (!identical(j$id, 1:length(j))) {
    stop("!identical(j$id, 1:length(j))")
  }
  rescue.id <- j$id[bRescue]

  if (!is.null(bam.ext.size)) {
    rescue.gr <- j[rescue.id, "id"]
    rescue.gr <- rescue.gr - ((width(rescue.gr)-1)/2)
    rescue.gr <- rescue.gr + bam.ext.size
  } else {
    rescue.gr <- sync.o$atac[rescue.id, ] %>% rowRanges() %>% mutate(id = rescue.id)
    stop("untested")
  }
  rescue.gr$id <- paste0("rescue_", rescue.gr$id)
  #### rescue function:
  samples.ori <- colnames(sync.o$atac)
  bam.map <- read.table(bam.sample.map, header = T)
  # attempt to rearrange:
  bam.map <- bam.map[match(samples.ori, bam.map$sample),]
  if (!identical(bam.map$sample, samples.ori)) {
    stop("sample.names must exactly match the colnames in atac.se. In the same order")
  }
  
  
  rescue.o <- bamFanc::cut.site.pipe(bam.files = bam.map$bam, region.gr = rescue.gr, name.col = "id", 
                                     sample.names = bam.map$sample, 
                                     ctrl.bed.files = bam.map$ctrl,
                                     shift = rescue.shift,
                                     report.dir = rescue.report.dir, tag = tag, 
                                     header.file = header.file, m.score.cutoff = m.score.cutoff, 
                                     threads.bam = threads, threads.region = 1)
  sync.o.rescue <- sync.o
  rescued.regions <- rescue.o$join$stat$region %>% sub("rescue_", "", .) %>% as.numeric()
  rescued.data <- rescue.o$join$stat %>% .[, colnames(.) != "region"] %>% as.matrix()
  if (!identical(colnames(sync.o.rescue$atac), colnames(rescued.data)))
    stop("colnames(sync.o.rescue$atac) != colnames(rescued.regions)")
  
  assays(sync.o.rescue$atac)$counts[rescued.regions,] <- rescued.data[,colnames(sync.o.rescue$atac)]
  ranges(rowRanges(sync.o.rescue$atac)[rescued.regions]) <- ranges(rescue.gr)
  # which((assays(sync.o.rescue$atac)$cpm[,1] - assays(sync.o$atac)$cpm[,1]) != 0)
  sync.o.rescue$atac <- se.fast.normalize(se = sync.o.rescue$atac, factor = factor, slot.name = "cpm")
  res <- list(sync.o = sync.o, sync.o.rescue = sync.o.rescue, rescue.o = rescue.o)
  return(res)
}

range.cor.2 <- function(sync.o, cage.samples = NULL, atac.samples = NULL, 
                        anno = "SYMBOL", transformation = function(x) return(log2(x + 1)), 
                        outfile = NULL, 
                        ...) {
  if (is.null(cage.samples))
    cage.samples <- colnames(sync.o$cage)
  if (is.null(atac.samples))
    atac.samples <- colnames(sync.o$atac)
  
  df <- data.frame(anno = mcols(sync.o$j)[, anno],
                   id = sync.o$j$id,
                   cage.id = sync.o$j$cage.id,
                   atac.id = sync.o$j$atac.id,
                   cage = assays(sync.o$cage)$tpm[, cage.samples] %>% apply(1, mean),
                   atac = assays(sync.o$atac)$cpm[, atac.samples] %>% apply(1, mean),
                   plotly = paste0(utilsFanc::gr.get.loci(sync.o$j), "..", anno))
  df <- df %>% mutate(anno.ez = gsub("[a-zA-Z]","",anno) %>% sub("\\-", "",.))
  df$anno.ez[!grepl("Klra\\d",df$anno)] <- NA
  
  p <- scFanc::xy.plot(df = df, x = "cage", y = "atac", plotly.var = "plotly", 
                       transformation = transformation,
                       highlight.var = "anno.ez", label.var = "anno.ez",
                       highlight.values = df$anno.ez %>% .[!is.na(.)], 
                       label.values = df$anno.ez %>% .[!is.na(.)], 
                       pt.size = 0.1, highlight.ptsize = 0.5, 
                       add.abline = F, nudge_x = 0.2, text.size = 3, 
                       outfile = outfile,
                       ...)
  
  return(list(p = p, df = df))
}
