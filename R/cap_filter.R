cap.filter.fanc <- function(gr, bam, 
                            thread.peak = 8, 
                            samtools = "/bar/cfan/anaconda2/envs/jupyter/bin/samtools") {
  
  
  if (!file.exist(paste0(bam, ".bai"))) {
    system(paste0(samtools, " index ", bam))
  }
  
  what <- c("flag" , "cigar", "seq")
  param <- Rsamtools::ScanBamParam(which = gr, what = what)
  reads.list <- Rsamtools::scanBam(bam, param = param, index = paste0(bam, ".bai"))
  reads.list <<- reads
  n.tss.sclip <- mclapply(reads.list, function(reads) {
    flags <- Rsamtools::bamFlagAsBitMatrix(reads$flag)
    
  }, mc.cores = thread.peak, mc.cleanup = T)
  
}

cap.filter.nakul <- function(df, ce, samples = NULL, plot.dir,
                             cap.pct.threads = 8) {
  system(paste0("mkdir -p ", plot.dir))
  if (is.null(samples))
    samples <- ce$sampleLabels
  else
    samples <- ce$sampleLabels %>% .[.%in% samples]

  # first create unannotated G ce.
  bams <- ce@colData %>% as.data.frame() %>% filter(sampleLabels %in% samples) %>% pull(inputFiles)
  G.CTSS <- paste0(bams, "_unannotatedG.CTSS")
  if (sum(!file.exists(G.CTSS)) > 0) {
    stop(paste0("these G.CTSS files do not exist: \n",
                paste0(G.CTSS[!file.exists(G.CTSS)], collapse = "\n")))
  }
  ce.G <- new("CAGEexp", genomeName = genomeName(ce), inputFiles = G.CTSS, inputFilesType = "ctss",
              sampleLabels = samples)
  getCTSS(ce.G)
  
  CTSSpeakcount <- lapply(list(ce, ce.G), function(ce.x) {
    ctss.tags <- CTSStagCount(object = ce.x)[, c("chr", "pos", "strand", samples)] %>%
      mutate(start = pos, end = pos) %>% mutate(pos = NULL) %>% 
      GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)  
    peakcount <- gr.get.coverage(gr.int = df, gr.signal = ctss.tags)
    return(peakcount)
  } )
  names(CTSSpeakcount) <- c("all", "uG")
  uG.pct <- CTSSpeakcount$uG/CTSSpeakcount$all 
  uG.pct[is.na(uG.pct)] <- 0
  colnames(uG.pct) <- paste0("cap.", colnames(uG.pct))
  
  ####################
  scFanc::rank.plot(df = uG.pct, vars = paste0("cap.",samples),
                    outfile = paste0(plot.dir, "/cap.pct.rank.png"))
  # ggplot(uG.pct, aes(x = cap.B6_Ly49Dn_mRNA)) + geom_density()
  qs <- lapply(samples, function(sample) {
    q <- quantile(uG.pct[, paste0("cap.", sample)], 0.1 * (1:10))
  }) %>% Reduce(cbind,.) %>% as.data.frame()
  colnames(qs) <- samples
  qs <- cbind(data.frame(quantile = rownames(qs)), qs)
  write.table(qs, paste0(plot.dir, "/quantile_table.tsv"), row.names = T, col.names = T, 
              sep = "\t", quote = F)
  ####################
  
  uG.pct$max <- apply(uG.pct, 1, max)
  uG.pct$min <- apply(uG.pct, 1, min)
  
  res <- list(uG.pct = uG.pct, CTSSpeakcount = CTSSpeakcount)
  return(res)
}

get.cap.daofeng <- function(bams, out.sh,
                            python2 = "/opt/apps/python2/bin/python2", 
                            py.script = "/bar/cfan/scripts/cage/bam2CTSS.py",
                            debug = F) {
  G.CTSS <- paste0(bams, "_unannotatedG.CTSS")
  cmds <- sapply(G.CTSS, function(x) {
    cmd <- paste0(python2, " ", py.script, " ", sub("_unannotatedG.CTSS$", "", x), " ", dirname(x))
    return(cmd)
  })
  cmds <- c("#!/bin/bash", cmds)
  system(paste0("mkdir -p ", dirname(out.sh)))
  write(cmds, out.sh, sep = "\n")
  system(paste0("chmod 755 ", out.sh))
  return(out.sh)
}

gr.get.coverage <- function(gr.int, gr.signal, samples = NULL) {
  if (is.null(samples))
    samples <- colnames(mcols(gr.signal))
  if (is.data.frame(gr.int))
    gr.int <- gr.int %>% GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = F)
  if (is.data.frame(gr.signal))
    gr.signal <- gr.signal %>% GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)
  # gr.int <- df %>% GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = F)
  # gr.signal <- ctss.tags %>% mutate(start = pos, end = pos) %>% mutate(pos = NULL) %>% 
  #   GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)
  
  o <- findOverlaps(query = gr.int, subject = gr.signal, ignore.strand = F) %>% as.data.frame()
  df.signal <- gr.signal %>% mcols() %>% as.data.frame()
  df.signal$subjectHits <- 1:nrow(df.signal)
  df.sum <- left_join(o, df.signal) %>% na.omit()
  df.sum <- df.sum %>% group_by(queryHits) %>% summarise_at(.vars = samples, .funs = sum) %>% 
    ungroup() %>% as.data.frame()
  df.sum <- df.sum %>% left_join(data.frame(queryHits = 1:length(gr.int)), .) 
  df.sum <- df.sum %>% mutate(queryHits = NULL)
  # t <- Reduce(rbind, resultapply)
  # t <- t %>% as.data.frame()
  # identical(as.numeric(t$B6_Ly49Dn_mRNA), as.numeric(df.sum$B6_Ly49Dn_mRNA))
  #returned TRUE. my fast function returns exactly the same result as Nakul's slow one
  return(df.sum)
}