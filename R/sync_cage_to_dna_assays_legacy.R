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

    # added on 2021-09-05 when re-running the scripts
    names(j) <- NULL
    #

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
    atac.se <- SummarizedExperiment::rbind(atac.se, re.read.se)
    # atac.se <- se.fast.normalize(se = atac.se, factor = factor, slot.name = "cpm")
  }

  sync.o <- list(j = j, cage = cage.se[j$cage.id, ],
                 atac = atac.se[j$atac.id, ])
  names(sync.o$j) <- NULL
  rownames(sync.o$cage) <- 1:nrow(sync.o$cage)
  rownames(sync.o$atac) <- 1:nrow(sync.o$atac)

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

range.cor.2 <- function(sync.o, df = NULL, cage.samples = NULL, atac.samples = NULL,
                        use.template.1 = F, use.template.2 = F,
                        anno = "SYMBOL", transformation = function(x) return(log2(x + 1)),
                        outfile = NULL, use.Ly49.names = T, return.df = F, add.cor = F,
                        ...) {
  if (is.null(df)) {
    if (is.null(cage.samples))
      cage.samples <- colnames(sync.o$cage)
    if (is.null(atac.samples))
      atac.samples <- colnames(sync.o$atac)
    names(sync.o$j) <- NULL
    rownames(sync.o$cage) <- 1:nrow(sync.o$cage)
    rownames(sync.o$atac) <- 1:nrow(sync.o$atac)
    df <- data.frame(anno = mcols(sync.o$j)[, anno],
                     id = sync.o$j$id,
                     cage.id = sync.o$j$cage.id,
                     atac.id = sync.o$j$atac.id,
                     cage = assays(sync.o$cage)$tpm[, cage.samples] %>% apply(1, mean),
                     atac = assays(sync.o$atac)$cpm[, atac.samples] %>% apply(1, mean),
                     plotly = paste0(utilsFanc::gr.get.loci(sync.o$j), "..", anno))

    df <- df %>% mutate(anno.ez = gsub("[a-zA-Z]","",anno) %>% sub("\\-", "",.))
    df$anno.ez[!grepl("Klra\\d",df$anno)] <- NA
    if (use.Ly49.names == T)
      df$anno.ez <- sapply(df$anno.ez, function(x) {
        if (!is.na(x))
          x <- letters[as.numeric(x)]
        return(x)
      })
    out.tsv <- tools::file_path_sans_ext(outfile) %>% paste0(".tsv")
    write.table(df, out.tsv, col.names = T, row.names = F, sep = "\t",
                quote = F)
  } else {
    if (is.character(df)) {
      df <- read.table(df, header = T)
    }
  }

  pro2.df <- df %>% filter(!is.na(anno.ez)) %>% filter(atac != max(atac)) %>%
    mutate(atac = transformation(atac), cage = transformation(cage))
  if (use.template.1) {
    p <- scFanc::xy.plot(df = df, x = "cage", y = "atac",
                         plotly.var = "anno", plotly.label.all = T,
                         transformation = transformation,
                         highlight.var = "anno.ez", label.var = "anno.ez",
                         highlight.values = df$anno.ez %>% .[!is.na(.)],
                         label.values = df$anno.ez %>% .[!is.na(.)],
                         pt.size = 0.1, highlight.ptsize = 0.12,
                         add.abline = F, nudge_x = 0.2, text.size = 3,
                         outfile = NULL, italic.label = T,
                         ...)
    if (add.cor == T) {
      p <- p + ggpubr::stat_cor(data = pro2.df, inherit.aes = F, mapping = aes(x = cage, y = atac))
      p <- p + geom_smooth(mapping = aes(x = cage, y = atac), data = pro2.df, linetype = "dashed", color = "green",
                           inherit.aes = F, method = "lm", se = F)
    }
    p <- scFanc::wrap.plots.fanc(plot.list = list(p), plot.out = outfile)
  } else if (use.template.2) {
    p <- scFanc::xy.plot(df = df, x = "cage", y = "atac",
                         plotly.var = "anno", plotly.label.all = T,
                         transformation = transformation,
                         highlight.var = "anno.ez", # label.var = "anno.ez",
                         highlight.values = df$anno.ez %>% .[!is.na(.)],
                         label.values = df$anno.ez %>% .[!is.na(.)],
                         pt.size = 0.0001, highlight.ptsize = 0.12,
                         add.abline = F, nudge_x = 0.2, text.size = 3,
                         outfile = NULL, italic.label = T,
                         ...)
    if (add.cor == T) {
      p <- p + ggpubr::stat_cor(
        data = pro2.df, inherit.aes = F,
        mapping = aes(x = cage, y = atac), size = 5 * 0.36,
        family = "Arial", cor.coef.name = "rho",
        label.x = 8.5, label.y = 0)
      p <- p + geom_smooth(
        mapping = aes(x = cage, y = atac), data = pro2.df,
        linetype = "dashed", color = "blue", size = 0.3,
        inherit.aes = F, method = "lm", se = F)
    }

    p <- p %>% utilsFanc::theme.fc.1(italic.x = F)

    ggsave(outfile, p, device = "pdf", width = 1.5, height = 1.5, dpi = 300)
    embedFonts(outfile)
  } else{
    p <- scFanc::xy.plot(df = df, x = "cage", y = "atac",
                         plotly.var = "anno", plotly.label.all = T,
                         transformation = transformation,
                         # highlight.var = "anno.ez", label.var = "anno.ez",
                         # highlight.values = df$anno.ez %>% .[!is.na(.)],
                         # label.values = df$anno.ez %>% .[!is.na(.)],
                         # pt.size = 0.1, highlight.ptsize = 0.12,
                         # add.abline = F, nudge_x = 0.2, text.size = 3,
                         # outfile = NULL, italic.label = T,
                         ...)
  }
  if (return.df == F)
    return(p)
  return(list(p = p, df = df))
}
