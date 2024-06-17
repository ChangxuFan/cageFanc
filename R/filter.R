strict.filter <- function(in.bams, general.cutoff, threads.master = 3, threads.sub = 4,
                          noMismatch = T, scfilter = T) {
  mclapply(in.bams, function(bam) {
    bam.noMismatch <- sub(".bam$", "_noMismatch.bam", bam)
    bam.scfilter <- sub(".bam$", "_scfilter.bam", bam.noMismatch)
    if (noMismatch == T) {
      stream.bam.core(bam.file = bam, fields = ALL.FIELDS, tags = CAGE.TAGS,header.file = NULL,
                      chunk.size = 30000, chunk.call.back = bam.filter.chunk,
                      chunk.call.back.param.list = list(outfile = bam.noMismatch,
                                                        threads = threads.sub,
                                                        filter.fun = bam.filter.no.mismatch,
                                                        filter.fun.params = list(tlen.cutoff = 10000)))
    }
    if (scfilter == T) {
      stream.bam.core(bam.file = bam.noMismatch, fields = ALL.FIELDS, tags = CAGE.TAGS,
                      chunk.size = 30000, chunk.call.back = bam.filter.chunk,
                      chunk.call.back.param.list = list(outfile = bam.scfilter,
                                                        # junkfile = "test/sth/mini_junk.bam",
                                                        # remove.mate = T,
                                                        filter.fun = bam.filter.softclip,
                                                        filter.fun.params = list(general.cutoff = general.cutoff,
                                                                                 threads = threads.sub)),
                      header.file = NULL)
    }

  }, mc.cores = threads.master, mc.cleanup = T)
}


umi.unique <- function(bams.in, out.dir = NULL, chunk.size = 100000, npar = 1) {
  print("this function uses Rsamtools::filterBam, which doesn't release memory when it's done")
  if (chunk.size > 4^9) {
    stop(paste0(
      "Your chunk size is probably too large. A 9 bp UMI can only label 4^9 (262144) transcripts"))
  }
  # note: by default, Rsamtools::filterBam process reads in 1M chunks.
  # this function is developed for each chunk. Each chunk doesn't see each other.
  # Not a big problem since if the files are alignment-position sorted.
  
  umi.unique.instance <- function(x) return(umi.unique.core(x = x))
  out.bams <- utilsFanc::safelapply(bams.in, function(bam) {
    if (is.null(out.dir)) {
      out.dir <- dirname(bam)
    }
    root.name <- basename(bam) %>% sub(".bam$", "", .)
    out.bam <- paste0(out.dir, "/", root.name, "_umiUnique.bam")
    Rsamtools::filterBam(file = Rsamtools::BamFile(bam, yieldSize = chunk.size), 
                         destination = out.bam,
                         filter = FilterRules(umi.unique.instance),
                         param = Rsamtools::ScanBamParam(tag = "RX", what = c("qname")))
    return(out.bam)
  }, threads = npar) %>% unlist()
  return(out.bams)
}

umi.unique.core <- function(x) {
  read.names <- x$qname[!duplicated(x$RX)]
  return(x$qname %in% read.names)
}
