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

star.fanc <- function(fastqs, genome.index, outdir, out.root.name="", thread, other.params.1 = ENCODE_RNA_STAR_PARAMS,
                      other.params.2="", stdout.file = NULL, run = T,
                      star.path = "/opt/apps/STAR/2.5.4b/STAR") {
  system(paste0("mkdir -p ", outdir))
  cmd <- paste0(star.path, " --genomeDir ", genome.index, " --readFilesIn ", paste0(fastqs, collapse = " "),
                " --runThreadN ", thread, " --outFileNamePrefix ", outdir, "/",out.root.name,
                " --outSAMtype BAM SortedByCoordinate")
  if (sum(grepl(".gz$", fastqs)) > 0) {
    cmd <- paste0(cmd, " --readFilesCommand zcat")
  }

  cmd <- paste0(cmd, " ", paste0(other.params.1, collapse = " "), " ", paste0(other.params.2, collapse = " "))
  utilsFanc::cmd.exec.fanc(cmd, stdout.file = stdout.file, intern = F, run = run)
  out.bam <- paste0(outdir, "/", out.root.name,"Aligned.sortedByCoord.out.bam")
  if (!file.exists(out.bam))
    stop(paste0(out.bam, " was not successfully generated"))
  return(out.bam)
}
