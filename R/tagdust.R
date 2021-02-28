xiaoyun.cage.block <- c("O:N", "S:ACACAG", "F:NNNNNNNN", "S:TATAGGG", "R:N")
juheon.cage.block <- c("O:N", "F:NNNNNNNNN", "S:TATAGGG", "R:N") # note: there are 9 Ns. Not 8.

tagdust <- function(fastqs, out.root.name, block.vec,
                    tagdust.path = "/bar/cfan/local/bin/bin/tagdust", thread,
                    show.finger.seq = T, other.params = "", run=T, stdout.file = NULL) {
  block.vec <- paste0("-", 1:length(block.vec), " ", block.vec)
  blocks <- paste0(block.vec, collapse = " ")
  cmd <- paste0(tagdust.path, " -t ", thread, " ", blocks, " ")
  if (show.finger.seq == T)
    cmd <- paste0(cmd, " -show_finger_seq")

  paste0(cmd, " ", other.params)
  cmd <- paste0(cmd, " ", paste0(fastqs, collapse = " "))
  cmd <- paste0(cmd, " -o ", out.root.name)

  utilsFanc::cmd.exec.fanc(cmd = cmd, stdout.file = stdout.file, intern = F, run = run)
  outfile <- paste0(out.root.name, "_READ",1:2,".fq")
  return(outfile)

}

tagdust.get.umi <- function(td.fastq, out.file = NULL, python2.path = "/opt/apps/python2/bin/python2",
                            run = T, stdout.file = NULL) {
  cat(paste0("tagdust.get.umi: a wrapper for python2 script ", "/bar/cfan/scripts/cage/TagdustOutput2UMI.py",
             ", which is copied from juheon's script ",
             "/scratch/jmaeng/TE_TF_cancer/TE_gene_chimericreads/script/TagdustOutput2UMI.py\n"))
  read1 <- td.fastq[grepl("READ1", td.fastq)]
  if (length(read1) != 1) {
    stop("something wrong with names of td.fastq")
  }

  if (is.null(out.file)) {
    out.file <- sub("fq", "UMI.fq", read1) %>% sub("fastq", "UMI.fastq", .)
  }
  cmd <- paste0(python2.path, " /bar/cfan/scripts/cage/TagdustOutput2UMI.py ", read1, " ", out.file)
  utilsFanc::cmd.exec.fanc(cmd = cmd, run = run, stdout.file = stdout.file, intern = F)
  if (!file.exists(out.file))
    stop(paste0(out.file, " failed to generate"))
  else
    return(out.file)
}
