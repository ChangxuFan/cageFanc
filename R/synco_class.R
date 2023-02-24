setClass("synco", 
         slots = c(j = "GRanges", 
                   cage.se = "SummarizedExperiment",
                   atac.se = "SummarizedExperiment"))

synco.gen.ez <- function(cage.se, atac.se, j = NULL, 
                         use.cage.as.j = T) {
  
  if (is.null(j)) {
    if (use.cage.as.j == T)
      j <- rowRanges(cage.se)
    else
      j <- rowRanges(atac.se)
  }
  
  if (length(j) != nrow(cage.se) | length(j) != nrow(atac.se)) {
    stop("for this function to work, the number of rows must be identical among" %>% 
           paste0(" j, cage.se and atac.se"))
  }
  
  rn <- names(j)
  if (is.null(rn)) {
    rn <- paste0("region_", 1:length(j))
  } else if (!any(is.na(suppressWarnings(as.numeric(rn))))) {
    rn <- paste0("region_", 1:length(j))
  }
  if (any(duplicated(rn))) {
    stop("rownames of j has duplicated elements")
  }
  
  names(j) <- rn
  rownames(cage.se) <- rn
  rownames(atac.se) <- rn
  synco <- new("synco", j = j, cage.se = cage.se, atac.se = atac.se)
  trash <- synco.validate(x = synco, raise.error = T)
  return(synco)
}


synco.validate <- function(x, raise.error = T, force.cpm = T) {
  if (!is(x, "synco")) {
    res <- list(bool = F, error= "not a synco")
  } else {
    rn <- names(x@j)
    scale.factor <- metadata(x@atac.se)$scale.factor
    if (length(rn) < 1) {
      res <- list(bool = F, error = "rownames of j is unset")
    } else if (any(duplicated(rn))) {
      res <- list(bool = F, error = "duplicated rownames exist in j")
    } else if (! identical(rn, rownames(x@atac.se)) | 
               ! identical(rn, rownames(x@cage.se))) {
      res <- list(bool = F, error = "rownames of j is not sync'ed to atac.se or cage.se")
    } else if (force.cpm == T && 
               any( abs(colSums(assays(x@atac.se)$cpm) - scale.factor) > 0.01 * scale.factor )) {
      res <- list(bool = F, error = "atac cpm does not conform to scale.factor in at least 1 sample")
    } else {
      res <- list(bool = T, error = NULL)
    }
  }
  if (raise.error == T && !is.null(res$error)) {
    stop(res$error)
  } 
  return(res)
}

synco.prep <- function(synco, annotate = T, rmbl = T, genome, 
                       genome.name = NULL, annoDb = NULL, TxDb = NULL) {
  trash <- synco.validate(x = synco, raise.error = T)
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
      annoDb <- "org.Hs.eg.db"
      TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    }
    anno <- synco@j %>%
      ChIPseeker::annotatePeak(TxDb=TxDb, annoDb = annoDb)
    if (length(anno@anno) != length(synco@j)) {
      stop ("length(anno@anno) != length(j)")
    }
    synco@j <- anno@anno
    rm(anno)
  }
  
  if (rmbl == T) {
    bl <- rtracklayer::import.bed(paste0("~/genomes/", genome, "/blacklist/", genome, ".blacklist.bed"))
    synco@j <- plyranges::filter_by_non_overlaps(x = synco@j, y = bl)
    synco <- synco.sync(synco)
  }
  trash <- synco.validate(x = synco, raise.error = T)
  return(synco)
}

"[.synco" <- function(x, i) {
  x@j <- x@j[i]
  x@atac.se <- x@atac.se[i,]
  x@cage.se <- x@cage.se[i,]
  x@atac.se <- x@atac.se %>% se.fast.normalize(factor = metadata(x@atac.se)$scale.factor, 
                                               slot = "cpm")
  trash <- synco.validate(x = x, raise.error = T)
  return(x)
}


synco.sync <- function(synco, rn = NULL) {
  if (is.null(rn)) {
    rn <- Reduce(intersect, list(names(synco@j), rownames(synco@atac.se), rownames(synco@cage.se)))
  }
  return(synco[rn])
}


synco.range.cor <- function(synco, cage.samples = NULL, atac.samples = NULL, use.template.2 = F,
                        anno = "SYMBOL", transformation = function(x) return(log2(x + 1)), 
                        outfile = NULL, use.Ly49.names = T, return.df = F, add.cor = F,
                        ...) {
  trash <- synco.validate(x = synco, raise.error = T, force.cpm = T)
  if (is.null(cage.samples))
    cage.samples <- colnames(synco@cage.se)
  if (is.null(atac.samples))
    atac.samples <- colnames(synco@atac.se)
  
  df <- data.frame(anno = mcols(synco@j)[, anno],
                   id = names(synco@j),
                   cage = assays(synco@cage.se)$tpm[, cage.samples, drop = F] %>% apply(1, mean),
                   atac = assays(synco@atac.se)$cpm[, atac.samples, drop = F] %>% apply(1, mean),
                   plotly = paste0(utilsFanc::gr.get.loci(synco@j), "..", anno))
  df <- df %>% mutate(anno.ez = gsub("[a-zA-Z]","",anno) %>% sub("\\-", "",.))
  df$anno.ez[!grepl("Klra\\d",df$anno)] <- NA
  if (use.Ly49.names == T)
    df$anno.ez <- sapply(df$anno.ez, function(x) {
      if (!is.na(x))
        x <- letters[as.numeric(x)]
      return(x)
    })
  pro2.df <- df %>% filter(!is.na(anno.ez)) %>% 
    mutate(atac = transformation(atac), cage = transformation(cage))
  
  if (use.template.2) {
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
  } else {
    p <- scFanc::xy.plot(df = df, x = "cage", y = "atac", plotly.var = "plotly", 
                         transformation = transformation,
                         highlight.var = "anno.ez", label.var = "anno.ez",
                         highlight.values = df$anno.ez %>% .[!is.na(.)], 
                         label.values = df$anno.ez %>% .[!is.na(.)], 
                         pt.size = 0.1, pt.color = "grey50", highlight.ptsize = 0.12, 
                         add.abline = F, nudge_x = 0.2, text.size = 3, 
                         outfile = outfile, italic.label = T,
                         ...)
    if (add.cor == T) {
      p <- p + ggpubr::stat_cor(data = pro2.df, inherit.aes = F, mapping = aes(x = cage, y = atac)) 
      p <- p + geom_smooth(mapping = aes(x = cage, y = atac), data = pro2.df, linetype = "dashed", color = "green",
                           inherit.aes = F, method = "lm", se = F) 
    }
    p <- scFanc::wrap.plots.fanc(plot.list = list(p), plot.out = outfile)
  }
  
  if (return.df == F)
    return(p)
  return(list(p = p, df = df))
}

synco.histone.pipe <- function(cage.se, atac.se = NULL, histone.samples.df, bSingleEnd,
                               buffer.up, buffer.down, genome,
                               cage.samples = NULL, atac.samples = NULL,
                               scale.factor = 10000,
                               root.name, out.dir) {
  if (!is.null(atac.se)) {
    if (is.character(atac.se))
      atac.se <- readRDS(atac.se)
  } else {
    atac.se <- create.histone.se(gr = rowRanges(cage.se), buffer.up = buffer.up, 
                                 buffer.down = buffer.down, samples.df = histone.samples.df, 
                                 bSingleEnd = bSingleEnd, scale.factor = scale.factor, 
                                 root.name = root.name, out.dir = out.dir)
  }

  synco <- synco.gen.ez(cage.se = cage.se, atac.se = atac.se, use.cage.as.j = T)
  synco <- synco.prep(synco = synco, annotate = T, rmbl = T, genome = genome)
  cor <- synco.range.cor(synco = synco, cage.samples = cage.samples, atac.samples = atac.samples, 
                         outfile = paste0(out.dir, "/", root.name, "_range_cor.png"))
  res <- list(synco = synco, cor = cor)
  return(res)
}

