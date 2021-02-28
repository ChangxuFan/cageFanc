peak.watching.quantile <- function(df, column, n.peaks, seed = 42, outfile = NULL) {
  if (is.character(df))
    df <- read.table(df, as.is = T, sep = "\t", header = T)
  df <- as.data.frame(df)
  set.seed(seed = seed)
  df.chosen <- df[sample(x = 1:nrow(df), size = n.peaks,replace = F),]
  df.chosen <- df.chosen %>% .[order(.[, column]),]
  if (!is.null(outfile))
    write.table(df.chosen, outfile, quote = F, col.names = T, row.names = F, sep = "\t")
  return(df.chosen)
}
