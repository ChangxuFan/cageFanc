ENCODE_RNA_STAR_PARAMS <- c(
  "--genomeLoad LoadAndKeep", "--outFilterMultimapNmax 20",
  "--alignSJoverhangMin 8", "--alignSJDBoverhangMin 1",
  "--outFilterMismatchNmax 999", "--outFilterMismatchNoverReadLmax 0.04",
  "--alignIntronMin 20", "--alignIntronMax 1000000", "--alignMatesGapMax 1000000",
  "--outSAMheaderHD @HD VN:1.4 SO:coordinate", "--outSAMunmapped Within",
  "--outFilterType BySJout", "--outSAMattributes NH HI AS NM MD",
  "--outSAMstrandField intronMotif",
  "--sjdbScore 1", "--limitBAMsortRAM 60000000000"
)

# what I left out from encode: "--outSAMheaderCommentFile COfile.txt", "--quantMode TranscriptomeSAM"
# mainly because I never use them.
# "--outSAMtype BAM SortedByCoordinate" was moved into the main function

star.fanc <- function (fastqs, genome.index, outdir, out.root.name = "", thread,
    other.params.1 = ENCODE_RNA_STAR_PARAMS, other.params.2 = "",
    stdout.file = NULL, run = T, star.path = STAR)
{
    system(paste0("mkdir -p ", outdir))
    cmd <- paste0(star.path, " --genomeDir ", genome.index, " --readFilesIn ",
        paste0(fastqs, collapse = " "), " --runThreadN ", thread,
        " --outFileNamePrefix ", outdir, "/", out.root.name,
        " --outSAMtype BAM SortedByCoordinate")
    if (sum(grepl(".gz$", fastqs)) > 0) {
        cmd <- paste0(cmd, " --readFilesCommand zcat")
    }
    cmd <- paste0(cmd, " ", paste0(other.params.1, collapse = " "),
        " ", paste0(other.params.2, collapse = " "))
    utilsFanc::cmd.exec.fanc(cmd, stdout.file = stdout.file,
        intern = F, run = run)
    out.bam <- paste0(outdir, "/", out.root.name, "Aligned.sortedByCoord.out.bam")
    if (!file.exists(out.bam))
        stop(paste0(out.bam, " was not successfully generated"))
    return(out.bam)
}
star.fanc.2 <- function(fastqs, suffix.regex = "_R*[12](_001)*.fastq.gz", genome.index, 
                      outdir, out.root.name=NULL,
                      other.params.1 = ENCODE_RNA_STAR_PARAMS,
                      other.params.2="", 
                      threads.sample = 1, threads.each = 1, 
                      stdout.file = NULL, run = T,
                      skip.align = F,
                      index = T, mkdup = T, flagstat = T, bw.gen = T, bw.norm = "RPKM",
                      filter.exp = "mapping_quality >= 30 and (not secondary_alignment) and (not duplicate)",
                      other.mkdup.options = " --hash-table-size=3000000 --overflow-list-size=3000000 ",
                      star.path = STAR,
                      samtools = SAMTOOLS) {
  
  # stop(paste0("just updated code to accept multiple samples. Haven't tested yet!!",
  #             "\nmain issue: out.root.name not working well because you need 1 for each sample", 
  #             "\nalso: indexing and bw has not been implemented"))
  system(paste0("mkdir -p ", outdir))
  roots <- fastqs %>% sub(suffix.regex, "", .)
  fastqs.list <- fastqs %>% split(., f = factor(roots, levels = unique(roots)))
  
  out.bams <-  utilsFanc::safelapply(names(fastqs.list), function(root) {
      fastq <- fastqs.list[[root]]
      root <- paste0(basename(root), "_")
      if (!is.null(out.root.name)) {
        if (length(fastq.list) != 1) {
          stop("you can only specify out.root.name when having exactly 1 sample")
        } else {
          root <- out.root.name
        }
      }
      if (!length(fastq) %in% c(1, 2)) {
        stop("!length(fastq) %in% c(1, 2)")
      }
      if (!skip.align) {
        cmd <- paste0(star.path, " --genomeDir ", genome.index, " --readFilesIn ", paste0(fastq, collapse = " "),
                " --runThreadN ", threads.each, " --outFileNamePrefix ", outdir, "/",root,
                " --outSAMtype BAM SortedByCoordinate")
        if (sum(grepl(".gz$", fastq)) > 0) {
          cmd <- paste0(cmd, " --readFilesCommand zcat")
        }

        cmd <- paste0(cmd, " ", paste0(other.params.1, collapse = " "), " ", paste0(other.params.2, collapse = " "))
        utilsFanc::cmd.exec.fanc(cmd, stdout.file = stdout.file, intern = F, run = run)
      }
      out.bam <- paste0(outdir, "/", root,"Aligned.sortedByCoord.out.bam")
      if (!file.exists(out.bam))
        stop(paste0(out.bam, " was not successfully generated"))
      if (index) {
        system(paste0(samtools, " index ", out.bam))
      }
      if (mkdup) {
        out.bam <- liteRnaSeqFanc::bam.dedup(in.bam = out.bam, thread = threads.each, remove = F, 
          other.options = other.mkdup.options)
      }
      if (flagstat) {
        flagstat.out <- paste0(tools::file_path_sans_ext(out.bam), "_flagstat.txt")
        system(paste0(samtools, " flagstat ", out.bam, " > ", flagstat.out))
      }
      if (!is.null(filter.exp)) {
        out.bam <- liteRnaSeqFanc::bam.filter(in.bam = out.bam, filter.expression = filter.exp,
                  thread = threads.each, creat.index = T)
      }
      if (bw.gen) {
        liteRnaSeqFanc::bam.browser(bam = out.bam, thread = threads.each, normalization = bw.norm)
      }
      return(out.bam)
    }, threads = threads.sample) %>% unlist()
  return(out.bams)
}
