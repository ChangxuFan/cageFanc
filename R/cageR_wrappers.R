bdg.bw.gen <- function(object, outdir, values) {
  # values should be raw or normalized
  wd <- getwd()
  system(paste0("mkdir -p ", outdir))
  setwd(outdir)
  try(exportCTSStoBedGraph(object, values = values, format = "BigWig"))
  try(exportCTSStoBedGraph(object, values = values, format = "bedGraph", oneFile = T))
  setwd(wd)
  return()
}

cager.cluster.2.granges <- function(df, keep.extra.columns = T) {
  return(GenomicRanges::makeGRangesFromDataFrame(df = df, keep.extra.columns = keep.extra.columns,
                                                 seqnames.field = "chr", start.field = "start", end.field = "end",
                                                 strand.field = "strand"))
}

