# this pipe is mainly adapted from Juheon's script:
#/scratch/jmaeng/TE_TF_cancer/TE_gene_chimericreads/script/pipe_align_nanoCAGE.sh

cage.pipe.juheon <- function(fastqs, root.name, out.dir="./", genome,
                             tagdust.block.vec = xiaoyun.cage.block,
                             cutadapt.min.read.length = 50,
                             thread, run=T) {
  if (length(fastqs) != 2)
    stop("length of fastqs must be 2")
  if (sum(!file.exists(fastqs))> 0)
    stop("both fastq files must exist")

  if (!grepl("/", genome))
    genome <- paste0("/bar/cfan/genomes/", genome, "/STAR_gencode_default")

  tagdust.dir <- paste0(out.dir, "/tagdust/",root.name, "/")
  bam.dir <- paste0(out.dir, "/star_aligned/", root.name, "/")

  system(paste0("mkdir -p ", tagdust.dir, " ", bam.dir))
  # note: this pipeline uses cagePipe as the prefix for intermediate variables. This is inheritated from test/mini.R.

  cagePipe.td <- tagdust(fastqs = fastqs, out.root.name = paste0(tagdust.dir, "/", root.name), block.vec = tagdust.block.vec,
                        thread = thread, show.finger.seq = T, run = run)
  cagePipe.umi <- tagdust.get.umi(td.fastq = cagePipe.td, run = run, stdout.file = paste0(tagdust.dir, "/get_umi.log"))

  cagePipe.td.trimmed <- cutadapt.fanc(fastqs = cagePipe.td, thread = thread, min.read.length = cutadapt.min.read.length,
                                      run = run, stdout.file = paste0(tagdust.dir, "/trim.log"))

  cagePipe.bam <- star.fanc(fastqs = cagePipe.td.trimmed, genome.index = genome, out.root.name = paste0(root.name, "_"),
                           outdir = bam.dir, thread = thread, run = run)
  cagePipe.bam.umi <- bam.add.umi(in.bam = cagePipe.bam, umi.fastq = cagePipe.umi, run = run, stdout.file = paste0(bam.dir, "/add_umi.log"))

  samtools.index(in.bam = cagePipe.bam.umi, thread = thread, run = run)

  cagePipe.bam.umi.dedup <- umi.tools.dedup.juheon(in.bam =  cagePipe.bam.umi, run = run, thread = thread)

  cagePipe.bam.umi.dedup.filtered <- samtools.filter.juheon(in.bam = cagePipe.bam.umi.dedup, thread = thread,run = run, index = T,
                                                           stdout.file = paste0(bam.dir, "/filter.log"))
  return(cagePipe.bam.umi.dedup.filtered)

}

para.ez <- function(sample.tsv, n.parallel=1, func, other.params) {
  if (is.character(sample.tsv))
    sample.tsv <- read.table(sample.tsv, as.is = T, header = T, sep = "\t", quote = "")

  res <- sample.tsv %>% utilsFanc::split.fanc("sample") %>%
    mclapply(function(x) {
      do.call(what = func, args = c(list(fastqs = x$fastq, root.name = x$sample[1]), other.params))
    }, mc.cores = n.parallel, mc.cleanup = T)
  return(res)
}





#
# t.f <- function(fastqs, root.name, a, b, c) {
#
#   return(paste0(fastqs, root.name, a, b ,c))
# }


# para.ez(data.frame(fastq = paste0("fastq", c(11,12,21,22)), sample = paste0("sample", c(1,1,2,2))),
#         n.parallel = 2, func = t.f, other.params = list(a = "mial", b="wang", c="huhuhu"))
#
#
# mclapply(X = c(1,2,3), FUN = t.f, b = "miao", c= "wang", mc.cores = 1)
