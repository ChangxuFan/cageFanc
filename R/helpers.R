fastq.mini.gen <- function(in.fastqs, nreads, out.fastqs, run = T) {
  if (length(in.fastqs) != length(out.fastqs))
    stop("the length of in.fastqs and out.fastqs must match")
  cat("note: gzip might print out 'stdout: Broken pipe', but it seems okay")
  trash <- mapply(function(x, y) {
    cmd <- paste0("cat", " ", x, " | head -", 4*nreads, " | gzip -nc  > ",y )
    if (grepl(".fastq.gz$", x)) {
      cmd <- paste0("z", cmd)
    }
    print(cmd)
    if (run == T) {
      system(paste0("mkdir -p ", dirname(y)))
      system(cmd)
    }

  }, in.fastqs, out.fastqs, SIMPLIFY = F)
  return()
}

bam.add.umi <- function(in.bam, umi.fastq, out.bam = NULL,
                        fgbio.path = "/opt/apps/java/jdk1.8.0_92/bin/java -jar /opt/apps/fgbio/1.0.0/fgbio-1.0.0.jar",
                        run = T, stdout.file = NULL) {
  if (is.null(out.bam))
    out.bam <- sub("bam$", "umi.bam", in.bam)
  if (in.bam == out.bam)
    stop("oops, out.bam turned out to be the same name as in.bam")
  cmd <- paste0(fgbio.path, " AnnotateBamWithUmis -i ", in.bam, " -f ", umi.fastq, " -o ", out.bam)
  utilsFanc::cmd.exec.fanc(cmd, intern = F, stdout.file = stdout.file, run = run)
  if (!file.exists(out.bam))
    stop(paste0(out.bam, " failed to generate"))
  return(out.bam)
}

samtools.index <- function(in.bam, thread, samtools.path = "/bar/cfan/anaconda2/envs/jupyter/bin/samtools",
                           run = T, stdout.file = NULL) {
  cmd <- paste0(samtools.path, " index -@ ", thread, " ",in.bam)
  utilsFanc::cmd.exec.fanc(cmd = cmd, stdout.file = stdout.file, run = run, intern = F)
  if (!file.exists(paste0(in.bam,".bai")))
    stop(paste0("the index for ", in.bam, " failed to generate"))
  return(in.bam)
}

umi.tools.dedup.juheon <- function(in.bam, out.bam=NULL, umi.tools.path = "/bar/cfan/anaconda2/envs/jupyter/bin/umi_tools",
                            run = T, index = T, stdout.file = NULL,
                            samtools.path = "/bar/cfan/anaconda2/envs/jupyter/bin/samtools",
                            thread = 1) {
  if (is.null(out.bam))
    out.bam <- sub("bam$", "dedup.bam", in.bam)
  if (in.bam == out.bam)
    stop("oops, out.bam turned out to be the same name as in.bam")
  out.root.name <- sub(".bam", "", out.bam)

  cmd <- paste0(umi.tools.path, " dedup --paired --extract-umi-method=tag --umi-tag=RX",
                " -I ", in.bam, " -L ", out.root.name, ".log",
                " -E ", out.root.name, ".error",
                " --output-stats=",out.root.name, ".stats",
                " -S ", out.bam)
  utilsFanc::cmd.exec.fanc(cmd, stdout.file = stdout.file, intern = F, run = run)

  if (!file.exists(out.bam))
    stop(paste0(out.bam, " failed to generate"))
  
  if (index == T) {
    samtools.index(in.bam = out.bam, thread = thread, run = run,
                   stdout.file = NULL, samtools.path = samtools.path)
  }
  
  return(out.bam)
}

samtools.filter.juheon <- function(in.bam, out.bam = NULL, thread,
                                   samtools.path = "/bar/cfan/anaconda2/envs/jupyter/bin/samtools",
                                   run = T, stdout.file = NULL, index = T) {
  if (is.null(out.bam))
    out.bam <- sub("bam$", "filtered.bam", in.bam)
  cmd <- paste0(samtools.path, " view -@ ", thread, " -f 0x2 -F 0x900 -q 255 ",
                " -bo ", out.bam, " ", in.bam)
  utilsFanc::cmd.exec.fanc(cmd, stdout.file = stdout.file, intern = F, run = run)

  if (!file.exists(out.bam))
    stop(paste0(out.bam, " failed to generate"))

  if (index == T) {
    samtools.index(in.bam = out.bam, thread = thread, run = run,
                   stdout.file = NULL, samtools.path = samtools.path)
  }
  return(out.bam)
}
